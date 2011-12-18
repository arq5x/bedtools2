// ***************************************************************************
// BamMultiReader_p.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Functionality for simultaneously reading multiple BAM files
// *************************************************************************

#include "api/BamAlignment.h"
#include "api/BamMultiReader.h"
#include "api/SamConstants.h"
#include "api/algorithms/Sort.h"
#include "api/internal/bam/BamMultiReader_p.h"
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
    : m_alignmentCache(0)
{ }

// dtor
BamMultiReaderPrivate::~BamMultiReaderPrivate(void) {
    Close();
}

// close all BAM files
bool BamMultiReaderPrivate::Close(void) {

    m_errorString.clear();

    if ( CloseFiles(Filenames()) )
        return true;
    else {
        const string currentError = m_errorString;
        const string message = string("error encountered while closing all files: \n\t") + currentError;
        SetErrorString("BamMultiReader::Close", message);
        return false;
    }
}

// close requested BAM file
bool BamMultiReaderPrivate::CloseFile(const string& filename) {

    m_errorString.clear();

    vector<string> filenames(1, filename);
    if ( CloseFiles(filenames) )
        return true;
    else {
        const string currentError = m_errorString;
        const string message = string("error while closing file: ") + filename + "\n" + currentError;
        SetErrorString("BamMultiReader::CloseFile", message);
        return false;
    }
}

// close requested BAM files
bool BamMultiReaderPrivate::CloseFiles(const vector<string>& filenames) {

    bool errorsEncountered = false;
    m_errorString.clear();

    // iterate over filenames
    vector<string>::const_iterator filesIter = filenames.begin();
    vector<string>::const_iterator filesEnd  = filenames.end();
    for ( ; filesIter != filesEnd; ++filesIter ) {
        const string& filename = (*filesIter);
        if ( filename.empty() ) continue;

        // iterate over readers
        vector<MergeItem>::iterator readerIter = m_readers.begin();
        vector<MergeItem>::iterator readerEnd  = m_readers.end();
        for ( ; readerIter != readerEnd; ++readerIter ) {
            MergeItem& item = (*readerIter);
            BamReader* reader = item.Reader;
            if ( reader == 0 ) continue;

            // if reader matches requested filename
            if ( reader->GetFilename() == filename ) {

                // remove reader's entry from alignment cache
                m_alignmentCache->Remove(reader);

                // clean up reader & its alignment
                if ( !reader->Close() ) {
                    m_errorString.append(1, '\t');
                    m_errorString.append(reader->GetErrorString());
                    m_errorString.append(1, '\n');
                    errorsEncountered = true;
                }
                delete reader;
                reader = 0;

                // delete reader's alignment entry
                BamAlignment* alignment = item.Alignment;
                delete alignment;
                alignment = 0;

                // remove reader from reader list
                m_readers.erase(readerIter);

                // on match, just go on to next filename
                // (no need to keep looking and item iterator is invalid now anyway)
                break;
            }
        }
    }

    // make sure alignment cache is cleaned up if all readers closed
    if ( m_readers.empty() && m_alignmentCache ) {
        m_alignmentCache->Clear();
        delete m_alignmentCache;
        m_alignmentCache = 0;
    }

    // return whether all readers closed OK
    return !errorsEncountered;
}

// creates index files for BAM files that don't have them
bool BamMultiReaderPrivate::CreateIndexes(const BamIndex::IndexType& type) {

    bool errorsEncountered = false;
    m_errorString.clear();

    // iterate over readers
    vector<MergeItem>::iterator itemIter = m_readers.begin();
    vector<MergeItem>::iterator itemEnd  = m_readers.end();
    for ( ; itemIter != itemEnd; ++itemIter ) {
        MergeItem& item = (*itemIter);
        BamReader* reader = item.Reader;
        if ( reader == 0 ) continue;

        // if reader doesn't have an index, create one
        if ( !reader->HasIndex() ) {
            if ( !reader->CreateIndex(type) ) {
                m_errorString.append(1, '\t');
                m_errorString.append(reader->GetErrorString());
                m_errorString.append(1, '\n');
                errorsEncountered = true;
            }
        }
    }

    // check for errors encountered before returning success/fail
    if ( errorsEncountered ) {
        const string currentError = m_errorString;
        const string message = string("error while creating index files: ") + "\n" + currentError;
        SetErrorString("BamMultiReader::CreateIndexes", message);
        return false;
    } else
        return true;
}

IMultiMerger* BamMultiReaderPrivate::CreateAlignmentCache(void) const {

    // fetch SamHeader
    SamHeader header = GetHeader();

    // if BAM files are sorted by position
    if ( header.SortOrder == Constants::SAM_HD_SORTORDER_COORDINATE )
        return new MultiMerger<Algorithms::Sort::ByPosition>();

    // if BAM files are sorted by read name
    if ( header.SortOrder == Constants::SAM_HD_SORTORDER_QUERYNAME )
        return new MultiMerger<Algorithms::Sort::ByName>();

    // otherwise "unknown" or "unsorted", use unsorted merger and just read in
    return new MultiMerger<Algorithms::Sort::Unsorted>();
}

const vector<string> BamMultiReaderPrivate::Filenames(void) const {

    // init filename container
    vector<string> filenames;
    filenames.reserve( m_readers.size() );

    // iterate over readers
    vector<MergeItem>::const_iterator itemIter = m_readers.begin();
    vector<MergeItem>::const_iterator itemEnd  = m_readers.end();
    for ( ; itemIter != itemEnd; ++itemIter ) {
        const MergeItem& item = (*itemIter);
        const BamReader* reader = item.Reader;
        if ( reader == 0 ) continue;

        // store filename if not empty
        const string& filename = reader->GetFilename();
        if ( !filename.empty() )
            filenames.push_back(filename);
    }

    // return result
    return filenames;
}

string BamMultiReaderPrivate::GetErrorString(void) const {
    return m_errorString;
}

SamHeader BamMultiReaderPrivate::GetHeader(void) const {
    const string& text = GetHeaderText();
    return SamHeader(text);
}

// makes a virtual, unified header for all the bam files in the multireader
string BamMultiReaderPrivate::GetHeaderText(void) const {

    // N.B. - right now, simply copies all header data from first BAM,
    //        and then appends RG's from other BAM files
    // TODO: make this more intelligent wrt other header lines/fields

    // if no readers open
    const size_t numReaders = m_readers.size();
    if ( numReaders == 0 ) return string();

    // retrieve first reader's header
    const MergeItem& firstItem = m_readers.front();
    const BamReader* reader = firstItem.Reader;
    if ( reader == 0 ) return string();
    SamHeader mergedHeader = reader->GetHeader();

    // iterate over any remaining readers (skipping the first)
    for ( size_t i = 1; i < numReaders; ++i ) {
        const MergeItem& item = m_readers.at(i);
        const BamReader* reader = item.Reader;
        if ( reader == 0 ) continue;

        // retrieve current reader's header
        const SamHeader currentHeader = reader->GetHeader();

        // append current reader's RG entries to merged header
        // N.B. - SamReadGroupDictionary handles duplicate-checking
        mergedHeader.ReadGroups.Add(currentHeader.ReadGroups);

        // TODO: merge anything else??
    }

    // return stringified header
    return mergedHeader.ToString();
}

// get next alignment among all files
bool BamMultiReaderPrivate::GetNextAlignment(BamAlignment& al) {
    return PopNextCachedAlignment(al, true);
}

// get next alignment among all files without parsing character data from alignments
bool BamMultiReaderPrivate::GetNextAlignmentCore(BamAlignment& al) {
    return PopNextCachedAlignment(al, false);
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
    if ( m_readers.empty() ) return 0;

    // return reference count from first reader
    const MergeItem& item = m_readers.front();
    const BamReader* reader = item.Reader;
    if ( reader == 0 ) return 0;
    else
        return reader->GetReferenceCount();
}

// returns vector of reference objects
const RefVector BamMultiReaderPrivate::GetReferenceData(void) const {

    // handle empty multireader
    if ( m_readers.empty() ) return RefVector();

    // return reference data from first BamReader
    const MergeItem& item = m_readers.front();
    const BamReader* reader = item.Reader;
    if ( reader == 0 ) return RefVector();
    else
        return reader->GetReferenceData();
}

// returns refID from reference name
int BamMultiReaderPrivate::GetReferenceID(const string& refName) const {

    // handle empty multireader
    if ( m_readers.empty() ) return -1;

    // return reference ID from first BamReader
    const MergeItem& item = m_readers.front();
    const BamReader* reader = item.Reader;
    if ( reader == 0 ) return -1;
    else
        return reader->GetReferenceID(refName);
}
// ---------------------------------------------------------------------------------------

// returns true if all readers have index data available
// this is useful to indicate whether Jump() or SetRegion() are possible
bool BamMultiReaderPrivate::HasIndexes(void) const {

    // handle empty multireader
    if ( m_readers.empty() )
        return false;

    bool result = true;

    // iterate over readers
    vector<MergeItem>::const_iterator readerIter = m_readers.begin();
    vector<MergeItem>::const_iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        const MergeItem& item = (*readerIter);
        const BamReader* reader = item.Reader;
        if ( reader  == 0 ) continue;

        // see if current reader has index data
        result &= reader->HasIndex();
    }

    return result;
}

// returns true if multireader has open readers
bool BamMultiReaderPrivate::HasOpenReaders(void) {

    // iterate over readers
    vector<MergeItem>::const_iterator readerIter = m_readers.begin();
    vector<MergeItem>::const_iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        const MergeItem& item = (*readerIter);
        const BamReader* reader = item.Reader;
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
    vector<MergeItem>::iterator readerIter = m_readers.begin();
    vector<MergeItem>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        MergeItem& item = (*readerIter);
        BamReader* reader = item.Reader;
        if ( reader == 0 ) continue;

        // jump in each BamReader to position of interest
        reader->Jump(refID, position);
    }

    // returns status of cache update
    return UpdateAlignmentCache();
}

// locate (& load) index files for BAM readers that don't already have one loaded
bool BamMultiReaderPrivate::LocateIndexes(const BamIndex::IndexType& preferredType) {

    bool errorsEncountered = false;
    m_errorString.clear();

    // iterate over readers
    vector<MergeItem>::iterator readerIter = m_readers.begin();
    vector<MergeItem>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        MergeItem& item = (*readerIter);
        BamReader* reader = item.Reader;
        if ( reader == 0 ) continue;

        // if reader has no index, try to locate one
        if ( !reader->HasIndex() ) {
            if ( !reader->LocateIndex(preferredType) ) {
                m_errorString.append(1, '\t');
                m_errorString.append(reader->GetErrorString());
                m_errorString.append(1, '\n');
                errorsEncountered = true;
            }
        }
    }

    // check for errors encountered before returning success/fail
    if ( errorsEncountered ) {
        const string currentError = m_errorString;
        const string message = string("error while locating index files: ") + "\n" + currentError;
        SetErrorString("BamMultiReader::LocatingIndexes", message);
        return false;
    } else
        return true;
}

// opens BAM files
bool BamMultiReaderPrivate::Open(const vector<string>& filenames) {

    m_errorString.clear();

    // put all current readers back at beginning (refreshes alignment cache)
    if ( !Rewind() ) {
        const string currentError = m_errorString;
        const string message = string("unable to rewind existing readers: \n\t") + currentError;
        SetErrorString("BamMultiReader::Open", message);
        return false;
    }

    // iterate over filenames
    bool errorsEncountered = false;
    vector<string>::const_iterator filenameIter = filenames.begin();
    vector<string>::const_iterator filenameEnd  = filenames.end();
    for ( ; filenameIter != filenameEnd; ++filenameIter ) {
        const string& filename = (*filenameIter);
        if ( filename.empty() ) continue;

        // attempt to open BamReader
        BamReader* reader = new BamReader;
        const bool readerOpened = reader->Open(filename);

        // if opened OK, store it
        if ( readerOpened )
            m_readers.push_back( MergeItem(reader, new BamAlignment) );

        // otherwise store error & clean up invalid reader
        else {
            m_errorString.append(1, '\t');
            m_errorString += string("unable to open file: ") + filename;
            m_errorString.append(1, '\n');
            errorsEncountered = true;

            delete reader;
            reader = 0;
        }
    }

    // check for errors while opening
    if ( errorsEncountered ) {
        const string currentError = m_errorString;
        const string message = string("unable to open all files: \t\n") + currentError;
        SetErrorString("BamMultiReader::Open", message);
        return false;
    }

    // check for BAM file consistency
    if ( !ValidateReaders() ) {
        const string currentError = m_errorString;
        const string message = string("unable to open inconsistent files: \t\n") + currentError;
        SetErrorString("BamMultiReader::Open", message);
        return false;
    }

    // update alignment cache
    return UpdateAlignmentCache();
}

bool BamMultiReaderPrivate::OpenFile(const std::string& filename) {
    vector<string> filenames(1, filename);
    if ( Open(filenames) )
        return true;
    else {
        const string currentError = m_errorString;
        const string message = string("could not open file: ") + filename + "\n\t" + currentError;
        SetErrorString("BamMultiReader::OpenFile", message);
        return false;
    }
}

bool BamMultiReaderPrivate::OpenIndexes(const vector<string>& indexFilenames) {

    // TODO: This needs to be cleaner - should not assume same order.
    //       And either way, shouldn't start at first reader.  Should start at
    //       first reader without an index?

    // make sure same number of index filenames as readers
    if ( m_readers.size() != indexFilenames.size() ) {
        const string message("size of index file list does not match current BAM file count");
        SetErrorString("BamMultiReader::OpenIndexes", message);
        return false;
    }

    bool errorsEncountered = false;
    m_errorString.clear();

    // iterate over BamReaders
    vector<string>::const_iterator indexFilenameIter = indexFilenames.begin();
    vector<string>::const_iterator indexFilenameEnd  = indexFilenames.end();
    vector<MergeItem>::iterator readerIter = m_readers.begin();
    vector<MergeItem>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        MergeItem& item = (*readerIter);
        BamReader* reader = item.Reader;

        // open index filename on reader
        if ( reader ) {
            const string& indexFilename = (*indexFilenameIter);
            if ( !reader->OpenIndex(indexFilename) ) {
                m_errorString.append(1, '\t');
                m_errorString += reader->GetErrorString();
                m_errorString.append(1, '\n');
                errorsEncountered = true;
            }
        }

        // increment filename iterator, skip if no more index files to open
        if ( ++indexFilenameIter == indexFilenameEnd )
            break;
    }

    // return success/fail
    if ( errorsEncountered ) {
        const string currentError = m_errorString;
        const string message = string("could not open all index files: \n\t") + currentError;
        SetErrorString("BamMultiReader::OpenIndexes", message);
        return false;
    } else
        return true;
}

bool BamMultiReaderPrivate::PopNextCachedAlignment(BamAlignment& al, const bool needCharData) {

    // skip if no alignments available
    if ( m_alignmentCache == 0 || m_alignmentCache->IsEmpty() )
        return false;

    // pop next merge item entry from cache
    MergeItem item = m_alignmentCache->TakeFirst();
    BamReader* reader = item.Reader;
    BamAlignment* alignment = item.Alignment;
    if ( reader == 0 || alignment == 0 )
        return false;

    // set char data if requested
    if ( needCharData ) {
        alignment->BuildCharData();
        alignment->Filename = reader->GetFilename();
    }

    // store cached alignment into destination parameter (by copy)
    al = *alignment;

    // load next alignment from reader & store in cache
    SaveNextAlignment(reader, alignment);
    return true;
}

// returns BAM file pointers to beginning of alignment data & resets alignment cache
bool BamMultiReaderPrivate::Rewind(void) {

    // skip if no readers open
    if ( m_readers.empty() )
        return true;

    // attempt to rewind files
    if ( !RewindReaders() ) {
        const string currentError = m_errorString;
        const string message = string("could not rewind readers: \n\t") + currentError;
        SetErrorString("BamMultiReader::Rewind", message);
        return false;
    }

    // return status of cache update
    return UpdateAlignmentCache();
}

// returns BAM file pointers to beginning of alignment data
bool BamMultiReaderPrivate::RewindReaders(void) {

    m_errorString.clear();
    bool errorsEncountered = false;

    // iterate over readers
    vector<MergeItem>::iterator readerIter = m_readers.begin();
    vector<MergeItem>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        MergeItem& item = (*readerIter);
        BamReader* reader = item.Reader;
        if ( reader == 0 ) continue;

        // attempt rewind on BamReader
        if ( !reader->Rewind() ) {
            m_errorString.append(1, '\t');
            m_errorString.append( reader->GetErrorString() );
            m_errorString.append(1, '\n');
            errorsEncountered = true;
        }
    }

    return !errorsEncountered;
}

void BamMultiReaderPrivate::SaveNextAlignment(BamReader* reader, BamAlignment* alignment) {

    // if can read alignment from reader, store in cache
    //
    // N.B. - lazy building of alignment's char data - populated only:
    //        automatically by alignment cache to maintain its sorting OR
    //        on demand from client call to future call to GetNextAlignment()

    if ( reader->GetNextAlignmentCore(*alignment) )
        m_alignmentCache->Add( MergeItem(reader, alignment) );
}

void BamMultiReaderPrivate::SetErrorString(const string& where, const string& what) const {
    static const string SEPARATOR = ": ";
    m_errorString = where + SEPARATOR + what;
}

bool BamMultiReaderPrivate::SetRegion(const BamRegion& region) {

    // NB: While it may make sense to track readers in which we can
    // successfully SetRegion, In practice a failure of SetRegion means "no
    // alignments here."  It makes sense to simply accept the failure,
    // UpdateAlignments(), and continue.

    // iterate over alignments
    vector<MergeItem>::iterator readerIter = m_readers.begin();
    vector<MergeItem>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        MergeItem& item = (*readerIter);
        BamReader* reader = item.Reader;
        if ( reader == 0 ) continue;

        // set region of interest
        reader->SetRegion(region);
    }

    // return status of cache update
    return UpdateAlignmentCache();
}

// updates our alignment cache
bool BamMultiReaderPrivate::UpdateAlignmentCache(void) {

    // create alignment cache if not created yet
    if ( m_alignmentCache == 0 ) {
        m_alignmentCache = CreateAlignmentCache();
        if ( m_alignmentCache == 0 ) {
            SetErrorString("BamMultiReader::UpdateAlignmentCache", "unable to create new alignment cache");
            return false;
        }
    }

    // clear any prior cache data
    m_alignmentCache->Clear();

    // iterate over readers
    vector<MergeItem>::iterator readerIter = m_readers.begin();
    vector<MergeItem>::iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        MergeItem& item = (*readerIter);
        BamReader* reader = item.Reader;
        BamAlignment* alignment = item.Alignment;
        if ( reader == 0 || alignment == 0 ) continue;

        // save next alignment from each reader in cache
        SaveNextAlignment(reader, alignment);
    }

    // if we get here, ok
    return true;
}

// ValidateReaders checks that all the readers point to BAM files representing
// alignments against the same set of reference sequences, and that the
// sequences are identically ordered.  If these checks fail the operation of
// the multireader is undefined, so we force program exit.
bool BamMultiReaderPrivate::ValidateReaders(void) const {

    m_errorString.clear();

    // skip if 0 or 1 readers opened
    if ( m_readers.empty() || (m_readers.size() == 1) )
        return true;

    // retrieve first reader
    const MergeItem& firstItem = m_readers.front();
    const BamReader* firstReader = firstItem.Reader;
    if ( firstReader == 0 ) return false;

    // retrieve first reader's header data
    const SamHeader& firstReaderHeader = firstReader->GetHeader();
    const string& firstReaderSortOrder = firstReaderHeader.SortOrder;

    // retrieve first reader's reference data
    const RefVector& firstReaderRefData = firstReader->GetReferenceData();
    const int firstReaderRefCount = firstReader->GetReferenceCount();
    const int firstReaderRefSize = firstReaderRefData.size();

    // iterate over all readers
    vector<MergeItem>::const_iterator readerIter = m_readers.begin();
    vector<MergeItem>::const_iterator readerEnd  = m_readers.end();
    for ( ; readerIter != readerEnd; ++readerIter ) {
        const MergeItem& item = (*readerIter);
        BamReader* reader = item.Reader;
        if ( reader == 0 ) continue;

        // get current reader's header data
        const SamHeader& currentReaderHeader = reader->GetHeader();
        const string& currentReaderSortOrder = currentReaderHeader.SortOrder;

        // check compatible sort order
        if ( currentReaderSortOrder != firstReaderSortOrder ) {
            const string message = string("mismatched sort order in ") + reader->GetFilename() +
                                   ", expected " + firstReaderSortOrder +
                                   ", but found " + currentReaderSortOrder;
            SetErrorString("BamMultiReader::ValidateReaders", message);
            return false;
        }

        // get current reader's reference data
        const RefVector currentReaderRefData = reader->GetReferenceData();
        const int currentReaderRefCount = reader->GetReferenceCount();
        const int currentReaderRefSize  = currentReaderRefData.size();

        // init reference data iterators
        RefVector::const_iterator firstRefIter   = firstReaderRefData.begin();
        RefVector::const_iterator firstRefEnd    = firstReaderRefData.end();
        RefVector::const_iterator currentRefIter = currentReaderRefData.begin();

        // compare reference counts from BamReader ( & container size, in case of BR error)
        if ( (currentReaderRefCount != firstReaderRefCount) ||
             (firstReaderRefSize    != currentReaderRefSize) )
        {
            stringstream s("");
            s << "mismatched reference count in " << reader->GetFilename()
              << ", expected " << firstReaderRefCount
              << ", but found " << currentReaderRefCount;
            SetErrorString("BamMultiReader::ValidateReaders", s.str());
            return false;
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
                stringstream s("");
                s << "mismatched references found in" << reader->GetFilename()
                  << "expected: " << endl;

                // print first reader's reference data
                RefVector::const_iterator refIter = firstReaderRefData.begin();
                RefVector::const_iterator refEnd  = firstReaderRefData.end();
                for ( ; refIter != refEnd; ++refIter ) {
                    const RefData& entry = (*refIter);
                    stringstream s("");
                    s << entry.RefName << " " << endl;
                }

                s << "but found: " << endl;

                // print current reader's reference data
                refIter = currentReaderRefData.begin();
                refEnd  = currentReaderRefData.end();
                for ( ; refIter != refEnd; ++refIter ) {
                    const RefData& entry = (*refIter);
                    s << entry.RefName << " " << entry.RefLength << endl;
                }

                SetErrorString("BamMultiReader::ValidateReaders", s.str());
                return false;
            }

            // update iterators
            ++firstRefIter;
            ++currentRefIter;
        }
    }

    // if we get here, everything checks out
    return true;
}
