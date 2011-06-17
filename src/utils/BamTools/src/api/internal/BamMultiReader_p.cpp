// ***************************************************************************
// BamMultiReader_p.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 April 2011 (DB)
// ---------------------------------------------------------------------------
// Functionality for simultaneously reading multiple BAM files
// *************************************************************************

#include <api/BamAlignment.h>
#include <api/BamMultiReader.h>
#include <api/internal/BamMultiMerger_p.h>
#include <api/internal/BamMultiReader_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
using namespace std;

// ctor
BamMultiReaderPrivate::BamMultiReaderPrivate(void)
    : m_alignments(0)
    , m_isCoreMode(false)
    , m_sortOrder(BamMultiReader::SortedByPosition)
{ }

// dtor
BamMultiReaderPrivate::~BamMultiReaderPrivate(void) {

    // close all open BAM readers
    Close();

    // clean up alignment cache
    delete m_alignments;
    m_alignments = 0;
}

// close all BAM files
void BamMultiReaderPrivate::Close(void) {
    CloseFiles( Filenames() );
}

// close requested BAM file
void BamMultiReaderPrivate::CloseFile(const string& filename) {    
    vector<string> filenames(1, filename);
    CloseFiles(filenames);
}

// close requested BAM files
void BamMultiReaderPrivate::CloseFiles(const vector<string>& filenames) {

    // iterate over filenames
    vector<string>::const_iterator filesIter = filenames.begin();
    vector<string>::const_iterator filesEnd  = filenames.end();
    for ( ; filesIter != filesEnd; ++filesIter ) {
        const string& filename = (*filesIter);
        if ( filename.empty() ) continue;

        // iterate over readers
        vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
        vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
        for ( ; readerIter != readerEnd; ++readerIter ) {
            BamReader* reader = (*readerIter).first;
            if ( reader == 0 ) continue;

            // if reader matches requested filename
            if ( reader->GetFilename() == filename ) {

                // remove reader/alignment from alignment cache
                m_alignments->Remove(reader);

                // close & delete reader
                reader->Close();
                delete reader;
                reader = 0;

                // delete reader's alignment entry
                BamAlignment* alignment = (*readerIter).second;
                delete alignment;
                alignment = 0;

                // remove reader from container
                m_readers.erase(readerIter);

                // on match, just go on to next filename
                // (no need to keep looking and iterator is invalid now anyway)
                break;
            }
        }
    }

    // make sure alignment cache is cleared if all readers are now closed
    if ( m_readers.empty() && m_alignments != 0 )
        m_alignments->Clear();
}

// creates index files for BAM files that don't have them
bool BamMultiReaderPrivate::CreateIndexes(const BamIndex::IndexType& type) {

    bool result = true;

    // iterate over readers
    vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // if reader doesn't have an index, create one
        if ( !reader->HasIndex() )
            result &= reader->CreateIndex(type);
    }

    return result;
}

IBamMultiMerger* BamMultiReaderPrivate::CreateMergerForCurrentSortOrder(void) const {
    switch ( m_sortOrder ) {
        case ( BamMultiReader::SortedByPosition ) : return new PositionMultiMerger;
        case ( BamMultiReader::SortedByReadName ) : return new ReadNameMultiMerger;
        case ( BamMultiReader::Unsorted )         : return new UnsortedMultiMerger;
        default :
            cerr << "BamMultiReader ERROR: requested sort order is unknown" << endl;
            return 0;
    }
}

const string BamMultiReaderPrivate::ExtractReadGroup(const string& headerLine) const {

    string readGroup("");
    stringstream headerLineSs(headerLine);
    string part;

    // parse @RG header line, looking for the ID: tag
    while( getline(headerLineSs, part, '\t') ) {
        stringstream partSs(part);
        string subtag;
        getline(partSs, subtag, ':');
        if ( subtag == "ID" ) {
            getline(partSs, readGroup, ':');
            break;
        }
    }
    return readGroup;
}

const vector<string> BamMultiReaderPrivate::Filenames(void) const {

    // init filename container
    vector<string> filenames;
    filenames.reserve( m_readers.size() );

    // iterate over readers
    vector<ReaderAlignment>::const_iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::const_iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        const BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // store filename if not empty
        const string filename = reader->GetFilename();
        if ( !filename.empty() )
            filenames.push_back( reader->GetFilename() );
    }

    // return result
    return filenames;
}

SamHeader BamMultiReaderPrivate::GetHeader(void) const {
    string text = GetHeaderText();
    return SamHeader(text);
}

// makes a virtual, unified header for all the bam files in the multireader
string BamMultiReaderPrivate::GetHeaderText(void) const {

    // TODO: merge SamHeader objects instead of parsing string data (again)

    // if only one reader is open
    if ( m_readers.size() == 1 ) {

        // just return reader's header text
        const ReaderAlignment& ra = m_readers.front();
        const BamReader* reader = ra.first;
        if ( reader ) return reader->GetHeaderText();

        // invalid reader
        return string();
    }

    string mergedHeader("");
    map<string, bool> readGroups;

    // foreach extraction entry (each BAM file)
    vector<ReaderAlignment>::const_iterator readerBegin = m_readers.begin();
    vector<ReaderAlignment>::const_iterator readerIter  = readerBegin;
    vector<ReaderAlignment>::const_iterator readerEnd   = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        const BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // get header from reader
        string headerText = reader->GetHeaderText();
        if ( headerText.empty() ) continue;

        // store header text in lines
        map<string, bool> currentFileReadGroups;
        const vector<string> lines = SplitHeaderText(headerText);

        // iterate over header lines
        vector<string>::const_iterator linesIter = lines.begin();
        vector<string>::const_iterator linesEnd  = lines.end();
        for ( ; linesIter != linesEnd; ++linesIter ) {

            // get next line from header, skip if empty
            const string headerLine = (*linesIter);
            if ( headerLine.empty() ) continue;

            // if first file, save HD & SQ entries
            // TODO: what if first file has empty header, should just check for empty 'mergedHeader' instead ?
            if ( readerIter == readerBegin ) {
                if ( headerLine.find("@HD") == 0 || headerLine.find("@SQ") == 0) {
                    mergedHeader.append(headerLine.c_str());
                    mergedHeader.append(1, '\n');
                }
            }

            // (for all files) append RG entries if they are unique
            if ( headerLine.find("@RG") == 0 ) {

                // extract read group name from line
                const string readGroup = ExtractReadGroup(headerLine);

                // make sure not to duplicate @RG entries
                if ( readGroups.find(readGroup) == readGroups.end() ) {
                    mergedHeader.append(headerLine.c_str() );
                    mergedHeader.append(1, '\n');
                    readGroups[readGroup] = true;
                    currentFileReadGroups[readGroup] = true;
                } else {
                    // warn iff we are reading one file and discover duplicated @RG tags in the header
                    // otherwise, we emit no warning, as we might be merging multiple BAM files with identical @RG tags
                    if ( currentFileReadGroups.find(readGroup) != currentFileReadGroups.end() ) {
                        cerr << "BamMultiReader WARNING: duplicate @RG tag " << readGroup
                             << " entry in header of " << reader->GetFilename() << endl;
                    }
                }
            }
        }
    }

    // return merged header text
    return mergedHeader;
}

// get next alignment among all files
bool BamMultiReaderPrivate::GetNextAlignment(BamAlignment& al) {
    m_isCoreMode = false;
    return LoadNextAlignment(al);
}

// get next alignment among all files without parsing character data from alignments
bool BamMultiReaderPrivate::GetNextAlignmentCore(BamAlignment& al) {
    m_isCoreMode = true;
    return LoadNextAlignment(al);
}

// ---------------------------------------------------------------------------------------
//
// NB: The following GetReferenceX() functions assume that we have identical
// references for all BAM files.  We enforce this by invoking the
// ValidateReaders() method to verify that our reference data is the same
// across all files on Open - so we will not encounter a situation in which
// there is a mismatch and we are still live.
//
// ---------------------------------------------------------------------------------------

// returns the number of reference sequences
int BamMultiReaderPrivate::GetReferenceCount(void) const {

    // handle empty multireader
    if ( m_readers.empty() )
        return 0;

    // return reference count from first reader
    const ReaderAlignment& ra = m_readers.front();
    const BamReader* reader = ra.first;
    if ( reader ) return reader->GetReferenceCount();

    // invalid reader
    return 0;
}

// returns vector of reference objects
const RefVector BamMultiReaderPrivate::GetReferenceData(void) const {

    // handle empty multireader
    if ( m_readers.empty() )
        return RefVector();

    // return reference data from first BamReader
    const ReaderAlignment& ra = m_readers.front();
    const BamReader* reader = ra.first;
    if ( reader ) return reader->GetReferenceData();

    // invalid reader
    return RefVector();
}

// returns refID from reference name
int BamMultiReaderPrivate::GetReferenceID(const string& refName) const {

    // handle empty multireader
    if ( m_readers.empty() )
        return -1;

    // return reference ID from first BamReader
    const ReaderAlignment& ra = m_readers.front();
    const BamReader* reader = ra.first;
    if ( reader ) return reader->GetReferenceID(refName);

    // invalid reader
    return -1;
}
// ---------------------------------------------------------------------------------------

// checks if any readers still have alignments
bool BamMultiReaderPrivate::HasAlignmentData(void) const {
    if ( m_alignments == 0 )
        return false;
    return !m_alignments->IsEmpty();
}

// returns true if all readers have index data available
// this is useful to indicate whether Jump() or SetRegion() are possible
bool BamMultiReaderPrivate::HasIndexes(void) const {

    // handle empty multireader
    if ( m_readers.empty() )
        return false;

    bool result = true;

    // iterate over readers
    vector<ReaderAlignment>::const_iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::const_iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        const BamReader* reader = (*readerIter).first;
        if ( reader  == 0 ) continue;

        // see if current reader has index data
        result &= reader->HasIndex();
    }

    return result;
}

// returns true if multireader has open readers
bool BamMultiReaderPrivate::HasOpenReaders(void) {

    // iterate over readers
    vector<ReaderAlignment>::const_iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::const_iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        const BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // return true whenever an open reader is found
        if ( reader->IsOpen() ) return true;
    }

    // no readers open
    return false;
}

// performs random-access jump using (refID, position) as a left-bound
bool BamMultiReaderPrivate::Jump(int refID, int position) {

    // NB: While it may make sense to track readers in which we can
    // successfully Jump, in practice a failure of Jump means "no
    // alignments here."  It makes sense to simply accept the failure,
    // UpdateAlignments(), and continue.

    // iterate over readers
    vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // attempt jump() on each
        if ( !reader->Jump(refID, position) ) {
            cerr << "BamMultiReader ERROR: could not jump " << reader->GetFilename()
                 << " to " << refID << ":" << position << endl;
        }
    }

    // update alignment cache & return success
    UpdateAlignmentCache();
    return true;
}

bool BamMultiReaderPrivate::LoadNextAlignment(BamAlignment& al) {

    // bail out if no more data to process
    if ( !HasAlignmentData() )
        return false;

    // "pop" next alignment and reader
    ReaderAlignment nextReaderAlignment = m_alignments->TakeFirst();
    BamReader* reader = nextReaderAlignment.first;
    BamAlignment* alignment = nextReaderAlignment.second;

    // store cached alignment into destination parameter (by copy)
    al = *alignment;

    // peek to next alignment & store in cache
    SaveNextAlignment(reader, alignment);

    // return success
    return true;
}

// locate (& load) index files for BAM readers that don't already have one loaded
bool BamMultiReaderPrivate::LocateIndexes(const BamIndex::IndexType& preferredType) {

    bool result = true;

    // iterate over readers
    vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // if reader has no index, try to locate one
        if ( !reader->HasIndex() )
            result &= reader->LocateIndex(preferredType);
    }

    return result;
}

// opens BAM files
bool BamMultiReaderPrivate::Open(const vector<string>& filenames) {

    // create alignment cache if neccessary
    if ( m_alignments == 0 ) {
        m_alignments = CreateMergerForCurrentSortOrder();
        if ( m_alignments == 0 ) return false;
    }

    // iterate over filenames
    vector<string>::const_iterator filenameIter = filenames.begin();
    vector<string>::const_iterator filenameEnd  = filenames.end();
    for ( ; filenameIter != filenameEnd; ++filenameIter ) {
        const string& filename = (*filenameIter);
        if ( filename.empty() ) continue;

        // attempt to open BamReader on filename
        BamReader* reader = OpenReader(filename);
        if ( reader == 0 ) continue;

        // store reader with new alignment
        m_readers.push_back( make_pair(reader, new BamAlignment) );
    }

    // validate & rewind any opened readers, also refreshes alignment cache
    if ( !m_readers.empty() ) {
        ValidateReaders();
        Rewind();
    }

    // return success
    return true;
}

bool BamMultiReaderPrivate::OpenFile(const std::string& filename) {
    vector<string> filenames(1, filename);
    return Open(filenames);
}

bool BamMultiReaderPrivate::OpenIndexes(const vector<string>& indexFilenames) {

    // TODO: This needs to be cleaner - should not assume same order.
    //       And either way, shouldn't start at first reader.  Should start at
    //       first reader without an index?

    // make sure same number of index filenames as readers
    if ( m_readers.size() != indexFilenames.size() || !indexFilenames.empty() )
        return false;

    // init result flag
    bool result = true;

    // iterate over BamReaders
    vector<string>::const_iterator indexFilenameIter = indexFilenames.begin();
    vector<string>::const_iterator indexFilenameEnd  = indexFilenames.end();
    vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        BamReader* reader = (*readerIter).first;

        // open index filename on reader
        if ( reader ) {
            const string& indexFilename = (*indexFilenameIter);
            result &= reader->OpenIndex(indexFilename);
        }

        // increment filename iterator, skip if no more index files to open
        if ( ++indexFilenameIter == indexFilenameEnd )
            break;
    }

    // TODO: validation ??

    // return success/fail
    return result;
}

BamReader* BamMultiReaderPrivate::OpenReader(const std::string& filename) {

    // create new BamReader
    BamReader* reader = new BamReader;

    // if reader opens OK
    if ( reader->Open(filename) ) {

        // attempt to read first alignment (sanity check)
        // if ok, then return BamReader pointer
        BamAlignment al;
        if ( reader->GetNextAlignmentCore(al) )
            return reader;

        // could not read alignment
        else {
            cerr << "BamMultiReader WARNING: Could not read first alignment from "
                 << filename << ", ignoring file" << endl;
        }
    }

    // reader could not open
    else {
        cerr << "BamMultiReader WARNING: Could not open "
              << filename << ", ignoring file" << endl;
    }

    // if we get here, there was a problem with this BAM file (opening or reading)
    // clean up memory allocation & return null pointer
    delete reader;
    return 0;
}

// print associated filenames to stdout
void BamMultiReaderPrivate::PrintFilenames(void) const {
    const vector<string>& filenames = Filenames();
    vector<string>::const_iterator filenameIter = filenames.begin();
    vector<string>::const_iterator filenameEnd  = filenames.end();
    for ( ; filenameIter != filenameEnd; ++filenameIter )
        cout << (*filenameIter) << endl;
}

// returns BAM file pointers to beginning of alignment data & resets alignment cache
bool BamMultiReaderPrivate::Rewind(void) {

    // clear out alignment cache
    m_alignments->Clear();

    // attempt to rewind files
    if ( !RewindReaders() ) {
        cerr << "BamMultiReader ERROR: could not rewind file(s) successfully";
        return false;
    }

    // reset cache & return success
    UpdateAlignmentCache();
    return true;
}

// returns BAM file pointers to beginning of alignment data
bool BamMultiReaderPrivate::RewindReaders(void) {

    bool result = true;

    // iterate over readers
    vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // attempt rewind on BamReader
        result &= reader->Rewind();
    }

    return result;
}

void BamMultiReaderPrivate::SaveNextAlignment(BamReader* reader, BamAlignment* alignment) {

    // must be in core mode && NOT sorting by read name to call GNACore()
    if ( m_isCoreMode && m_sortOrder != BamMultiReader::SortedByReadName ) {
        if ( reader->GetNextAlignmentCore(*alignment) )
            m_alignments->Add( make_pair(reader, alignment) );
    }

    // not in core mode and/or sorting by readname, must call GNA()
    else {
        if ( reader->GetNextAlignment(*alignment) )
            m_alignments->Add( make_pair(reader, alignment) );
    }
}

// sets the index caching mode on the readers
void BamMultiReaderPrivate::SetIndexCacheMode(const BamIndex::IndexCacheMode mode) {

    // iterate over readers
    vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // set reader's index cache mode
        reader->SetIndexCacheMode(mode);
    }
}

bool BamMultiReaderPrivate::SetRegion(const BamRegion& region) {

    // NB: While it may make sense to track readers in which we can
    // successfully SetRegion, In practice a failure of SetRegion means "no
    // alignments here."  It makes sense to simply accept the failure,
    // UpdateAlignments(), and continue.

    // iterate over alignments
    vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;

        // attempt to set BamReader's region of interest
        if ( !reader->SetRegion(region) ) {
            cerr << "BamMultiReader WARNING: could not jump " << reader->GetFilename() << " to "
                 << region.LeftRefID  << ":" << region.LeftPosition   << ".."
                 << region.RightRefID << ":" << region.RightPosition  << endl;
        }
    }

    // update alignment cache & return success
    UpdateAlignmentCache();
    return true;
}

void BamMultiReaderPrivate::SetSortOrder(const BamMultiReader::SortOrder& order) {

    // skip if no change needed
    if ( m_sortOrder == order ) return;

    // set new sort order
    m_sortOrder = order;

    // create new alignment cache based on sort order
    IBamMultiMerger* newAlignmentCache = CreateMergerForCurrentSortOrder();
    if ( newAlignmentCache == 0 ) return; // print error?

    // copy old cache contents to new cache
    while ( m_alignments->Size() > 0 ) {
        ReaderAlignment value = m_alignments->TakeFirst(); // retrieves & 'pops'
        newAlignmentCache->Add(value);
    }

    // remove old cache structure & point to new cache
    delete m_alignments;
    m_alignments = newAlignmentCache;
}

// splits the entire header into a list of strings
const vector<string> BamMultiReaderPrivate::SplitHeaderText(const string& headerText) const {

    stringstream header(headerText);
    string item;

    vector<string> lines;
    while ( getline(header, item) )
        lines.push_back(item);
    return lines;
}

// updates our alignment cache
void BamMultiReaderPrivate::UpdateAlignmentCache(void) {

    // skip if invalid alignment cache
    if ( m_alignments == 0 ) return;

    // clear the cache
    m_alignments->Clear();

    // seed cache with fully-populated alignments
    // further updates will fill with full/core-only as requested
    m_isCoreMode = false;

    // iterate over readers
    vector<ReaderAlignment>::iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        BamReader* reader = (*readerIter).first;
        BamAlignment* alignment = (*readerIter).second;
        if ( reader == 0 || alignment == 0 ) continue;

        // save next alignment from each reader in cache
        SaveNextAlignment(reader, alignment);
    }
}

// ValidateReaders checks that all the readers point to BAM files representing
// alignments against the same set of reference sequences, and that the
// sequences are identically ordered.  If these checks fail the operation of
// the multireader is undefined, so we force program exit.
void BamMultiReaderPrivate::ValidateReaders(void) const {

    // retrieve first reader data
    const BamReader* firstReader = m_readers.front().first;
    if ( firstReader == 0 ) return;
    const RefVector firstReaderRefData = firstReader->GetReferenceData();
    const int firstReaderRefCount = firstReader->GetReferenceCount();
    const int firstReaderRefSize = firstReaderRefData.size();

    // iterate over all readers
    vector<ReaderAlignment>::const_iterator readerIter = m_readers.begin();
    vector<ReaderAlignment>::const_iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {

        // get current reader data
        BamReader* reader = (*readerIter).first;
        if ( reader == 0 ) continue;
        const RefVector currentReaderRefData = reader->GetReferenceData();
        const int currentReaderRefCount = reader->GetReferenceCount();
        const int currentReaderRefSize  = currentReaderRefData.size();

        // init container iterators
        RefVector::const_iterator firstRefIter   = firstReaderRefData.begin();
        RefVector::const_iterator firstRefEnd    = firstReaderRefData.end();
        RefVector::const_iterator currentRefIter = currentReaderRefData.begin();

        // compare reference counts from BamReader ( & container size, in case of BR error)
        if ( (currentReaderRefCount != firstReaderRefCount) ||
             (firstReaderRefSize    != currentReaderRefSize) )
        {
            cerr << "BamMultiReader ERROR: mismatched number of references in " << reader->GetFilename()
                 << " expected " << firstReaderRefCount
                 << " reference sequences but only found " << currentReaderRefCount << endl;
            exit(1);
        }

        // this will be ok; we just checked above that we have identically-sized sets of references
        // here we simply check if they are all, in fact, equal in content
        while ( firstRefIter != firstRefEnd ) {
            const RefData& firstRef   = (*firstRefIter);
            const RefData& currentRef = (*currentRefIter);

            // compare reference name & length
            if ( (firstRef.RefName   != currentRef.RefName) ||
                 (firstRef.RefLength != currentRef.RefLength) )
            {
                cerr << "BamMultiReader ERROR: mismatched references found in " << reader->GetFilename()
                     << " expected: " << endl;

                // print first reader's reference data
                RefVector::const_iterator refIter = firstReaderRefData.begin();
                RefVector::const_iterator refEnd  = firstReaderRefData.end();
                for ( ; refIter != refEnd; ++refIter ) {
                    const RefData& entry = (*refIter);
                    cerr << entry.RefName << " " << entry.RefLength << endl;
                }

                cerr << "but found: " << endl;

                // print current reader's reference data
                refIter = currentReaderRefData.begin();
                refEnd  = currentReaderRefData.end();
                for ( ; refIter != refEnd; ++refIter ) {
                    const RefData& entry = (*refIter);
                    cerr << entry.RefName << " " << entry.RefLength << endl;
                }

                exit(1);
            }

            // update iterators
            ++firstRefIter;
            ++currentRefIter;
        }
    }
}
