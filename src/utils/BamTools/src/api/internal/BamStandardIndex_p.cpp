// ***************************************************************************
// BamStandardIndex.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 16 June 2011 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the standardized BAM index format (".bai")
// ***************************************************************************

#include <api/BamAlignment.h>
#include <api/internal/BamReader_p.h>
#include <api/internal/BamStandardIndex_p.h>
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>
using namespace std;

// static BamStandardIndex constants
const int BamStandardIndex::MAX_BIN               = 37450;  // =(8^6-1)/7+1
const int BamStandardIndex::BAM_LIDX_SHIFT        = 14;
const string BamStandardIndex::BAI_EXTENSION      = ".bai";
const char* const BamStandardIndex::BAI_MAGIC     = "BAI\1";
const int BamStandardIndex::SIZEOF_ALIGNMENTCHUNK = sizeof(uint64_t)*2;
const int BamStandardIndex::SIZEOF_BINCORE        = sizeof(uint32_t) + sizeof(int32_t);
const int BamStandardIndex::SIZEOF_LINEAROFFSET   = sizeof(uint64_t);

// ctor
BamStandardIndex::BamStandardIndex(Internal::BamReaderPrivate* reader)
    : BamIndex(reader)
    , m_indexStream(0)
    , m_cacheMode(BamIndex::LimitedIndexCaching)
    , m_buffer(0)
    , m_bufferLength(0)
{
     m_isBigEndian = BamTools::SystemIsBigEndian();
}

// dtor
BamStandardIndex::~BamStandardIndex(void) {
    CloseFile();
}

bool BamStandardIndex::AdjustRegion(const BamRegion& region, uint32_t& begin, uint32_t& end) {

    // retrieve references from reader
    const RefVector& references = m_reader->GetReferenceData();

    // make sure left-bound position is valid
    if ( region.LeftPosition > references.at(region.LeftRefID).RefLength )
        return false;

    // set region 'begin'
    begin = (unsigned int)region.LeftPosition;

    // if right bound specified AND left&right bounds are on same reference
    // OK to use right bound position as region 'end'
    if ( region.isRightBoundSpecified() && ( region.LeftRefID == region.RightRefID ) )
        end = (unsigned int)region.RightPosition;

    // otherwise, set region 'end' to last reference base
    else end = (unsigned int)references.at(region.LeftRefID).RefLength - 1;

    // return success
    return true;
}

void BamStandardIndex::CalculateCandidateBins(const uint32_t& begin,
                                              const uint32_t& end,
                                              set<uint16_t>& candidateBins)
{
    // initialize list, bin '0' is always a valid bin
    candidateBins.insert(0);

    // get rest of bins that contain this region
    unsigned int k;
    for (k =    1 + (begin>>26); k <=    1 + (end>>26); ++k) { candidateBins.insert(k); }
    for (k =    9 + (begin>>23); k <=    9 + (end>>23); ++k) { candidateBins.insert(k); }
    for (k =   73 + (begin>>20); k <=   73 + (end>>20); ++k) { candidateBins.insert(k); }
    for (k =  585 + (begin>>17); k <=  585 + (end>>17); ++k) { candidateBins.insert(k); }
    for (k = 4681 + (begin>>14); k <= 4681 + (end>>14); ++k) { candidateBins.insert(k); }
}

bool BamStandardIndex::CalculateCandidateOffsets(const BaiReferenceSummary& refSummary,
                                                 const uint64_t& minOffset,
                                                 set<uint16_t>& candidateBins,
                                                 vector<int64_t>& offsets)
{
    // attempt seek to first bin
    if ( !Seek(refSummary.FirstBinFilePosition, SEEK_SET) )
        return false;

    // iterate over reference bins
    uint32_t binId;
    int32_t numAlignmentChunks;
    set<uint16_t>::iterator candidateBinIter;
    for ( int i = 0; i < refSummary.NumBins; ++i ) {

        // read bin contents (if successful, alignment chunks are now in m_buffer)
        if ( !ReadBinIntoBuffer(binId, numAlignmentChunks) )
            return false;

        // see if bin is a 'candidate bin'
        candidateBinIter = candidateBins.find(binId);

        // if not, move on to next bin
        if ( candidateBinIter == candidateBins.end() )
            continue;

        // otherwise, check bin's contents against for overlap
        else {

            unsigned int offset = 0;
            uint64_t chunkStart;
            uint64_t chunkStop;

            // iterate over alignment chunks
            for (int j = 0; j < numAlignmentChunks; ++j ) {

                // read chunk start & stop from buffer
                memcpy((char*)&chunkStart, m_buffer+offset, sizeof(uint64_t));
                offset += sizeof(uint64_t);
                memcpy((char*)&chunkStop, m_buffer, sizeof(uint64_t));
                offset += sizeof(uint64_t);

                // swap endian-ness if necessary
                if ( m_isBigEndian ) {
                    SwapEndian_64(chunkStart);
                    SwapEndian_64(chunkStop);
                }

                // store alignment chunk's start offset
                // if its stop offset is larger than our 'minOffset'
                if ( chunkStop >= minOffset )
                    offsets.push_back(chunkStart);
            }

            // 'pop' bin ID from candidate bins set
            candidateBins.erase(candidateBinIter);

            // quit if no more candidates
            if ( candidateBins.empty() )
                break;
        }
    }

    // return success
    return true;
}

uint64_t BamStandardIndex::CalculateMinOffset(const BaiReferenceSummary& refSummary,
                                              const uint32_t& begin)
{
    // if no linear offsets exist, return 0
    if ( refSummary.NumLinearOffsets == 0 )
        return 0;

    // if 'begin' starts beyond last linear offset, use the last linear offset as minimum
    // else use the offset corresponding to the requested start position
    const int shiftedBegin = begin>>BamStandardIndex::BAM_LIDX_SHIFT;
    if ( shiftedBegin >= refSummary.NumLinearOffsets )
        return LookupLinearOffset( refSummary, refSummary.NumLinearOffsets-1 );
    else
        return LookupLinearOffset( refSummary, shiftedBegin );
}

void BamStandardIndex::CheckBufferSize(char*& buffer,
                                       unsigned int& bufferLength,
                                       const unsigned int& requestedBytes)
{
    try {
        if ( requestedBytes > bufferLength ) {
            bufferLength = requestedBytes + 10;
            delete[] buffer;
            buffer = new char[bufferLength];
        }
    } catch ( std::bad_alloc ) {
        cerr << "BamStandardIndex ERROR: out of memory when allocating "
             << requestedBytes << " byes" << endl;
        exit(1);
    }
}

void BamStandardIndex::CheckBufferSize(unsigned char*& buffer,
                                       unsigned int& bufferLength,
                                       const unsigned int& requestedBytes)
{
    try {
        if ( requestedBytes > bufferLength ) {
            bufferLength = requestedBytes + 10;
            delete[] buffer;
            buffer = new unsigned char[bufferLength];
        }
    } catch ( std::bad_alloc ) {
        cerr << "BamStandardIndex ERROR: out of memory when allocating "
             << requestedBytes << " byes" << endl;
        exit(1);
    }
}

bool BamStandardIndex::CheckMagicNumber(void) {

    // check 'magic number' to see if file is BAI index
    char magic[4];
    size_t elementsRead = fread(magic, sizeof(char), 4, m_indexStream);
    if ( elementsRead != 4 ) {
        cerr << "BamStandardIndex ERROR: could not read format 'magic number'" << endl;
        return false;
    }

    // compare to expected value
    if ( strncmp(magic, BamStandardIndex::BAI_MAGIC, 4) != 0 ) {
        cerr << "BamStandardIndex ERROR: invalid format" << endl;
        return false;
    }

    // otherwise OK
    return true;
}

void BamStandardIndex::ClearReferenceEntry(BaiReferenceEntry& refEntry) {
    refEntry.ID = -1;
    refEntry.Bins.clear();
    refEntry.LinearOffsets.clear();
}

void BamStandardIndex::CloseFile(void) {

    // close file stream
    if ( IsFileOpen() )
        fclose(m_indexStream);

    // clear index file summary data
    m_indexFileSummary.clear();

    // clean up I/O buffer
    delete[] m_buffer;
    m_buffer = 0;
    m_bufferLength = 0;
}

// builds index from associated BAM file & writes out to index file
bool BamStandardIndex::Create(void) {

    // return false if BamReader is invalid or not open
    if ( m_reader == 0 || !m_reader->IsOpen() ) {
        cerr << "BamStandardIndex ERROR: BamReader is not open"
             << ", aborting index creation" << endl;
        return false;
    }

    // rewind BamReader
    if ( !m_reader->Rewind() ) {
        cerr << "BamStandardIndex ERROR: could not rewind BamReader to create index"
             << ", aborting index creation" << endl;
        return false;
    }

    // open new index file (read & write)
    string indexFilename = m_reader->Filename() + Extension();
    if ( !OpenFile(indexFilename, "w+b") ) {
        cerr << "BamStandardIndex ERROR: could not open ouput index file: " << indexFilename
             << ", aborting index creation" << endl;
        return false;
    }

    // initialize BaiFileSummary with number of references
    const int& numReferences = m_reader->GetReferenceCount();
    ReserveForSummary(numReferences);

    // initialize output file
    bool createdOk = true;
    createdOk &= WriteHeader();

    // set up bin, ID, offset, & coordinate markers
    const uint32_t defaultValue = 0xffffffffu;
    uint32_t currentBin    = defaultValue;
    uint32_t lastBin       = defaultValue;
    int32_t  currentRefID  = defaultValue;
    int32_t  lastRefID     = defaultValue;
    uint64_t currentOffset = (uint64_t)m_reader->Tell();
    uint64_t lastOffset    = currentOffset;
    int32_t  lastPosition  = defaultValue;

    // iterate through alignments in BAM file
    BamAlignment al;
    BaiReferenceEntry refEntry;
    while ( m_reader->LoadNextAlignment(al) ) {

        // changed to new reference
        if ( lastRefID != al.RefID ) {

            // if not first reference, save previous reference data
            if ( lastRefID != (int32_t)defaultValue ) {

                SaveAlignmentChunkToBin(refEntry.Bins, currentBin, currentOffset, lastOffset);
                createdOk &= WriteReferenceEntry(refEntry);
                ClearReferenceEntry(refEntry);

                // write any empty references between (but *NOT* including) lastRefID & al.RefID
                for ( int i = lastRefID+1; i < al.RefID; ++i ) {
                    BaiReferenceEntry emptyEntry(i);
                    createdOk &= WriteReferenceEntry(emptyEntry);
                }

                // update bin markers
                currentOffset = lastOffset;
                currentBin    = al.Bin;
                lastBin       = al.Bin;
                currentRefID  = al.RefID;
            }

            // first pass
            // write any empty references up to (but *NOT* including) al.RefID
            else {
                for ( int i = 0; i < al.RefID; ++i ) {
                    BaiReferenceEntry emptyEntry(i);
                    createdOk &= WriteReferenceEntry(emptyEntry);
                }
            }

            // update reference markers
            refEntry.ID = al.RefID;
            lastRefID   = al.RefID;
            lastBin     = defaultValue;
        }

        // if lastPosition greater than current alignment position - file not sorted properly
        else if ( lastPosition > al.Position ) {
            cerr << "BamStandardIndex ERROR: BAM file is not properly sorted by coordinate"
                 << ", aborting index creation"
                 << endl
                 << "At alignment: " << al.Name
                 << " : previous position " << lastPosition
                 << " > this alignment position " << al.Position
                 << " on reference id: " << al.RefID << endl;
            return false;
        }

        // if alignment's ref ID is valid & its bin is not a 'leaf'
        if ( (al.RefID >= 0) && (al.Bin < 4681) )
            SaveLinearOffsetEntry(refEntry.LinearOffsets, al.Position, al.GetEndPosition(), lastOffset);

        // changed to new BAI bin
        if ( al.Bin != lastBin ) {

            // if not first bin on reference, save previous bin data
            if ( currentBin != defaultValue )
                SaveAlignmentChunkToBin(refEntry.Bins, currentBin, currentOffset, lastOffset);

            // update markers
            currentOffset = lastOffset;
            currentBin    = al.Bin;
            lastBin       = al.Bin;
            currentRefID  = al.RefID;

            // if invalid RefID, break out
            if ( currentRefID < 0 )
                break;
        }

        // make sure that current file pointer is beyond lastOffset
        if ( m_reader->Tell() <= (int64_t)lastOffset ) {
            cerr << "BamStandardIndex ERROR: calculating offsets failed"
                 << ", aborting index creation" << endl;
            return false;
        }

        // update lastOffset & lastPosition
        lastOffset   = m_reader->Tell();
        lastPosition = al.Position;
    }

    // after finishing alignments, if any data was read, check:
    if ( currentRefID >= 0 ) {

        // store last alignment chunk to its bin, then write last reference entry with data
        SaveAlignmentChunkToBin(refEntry.Bins, currentBin, currentOffset, lastOffset);
        createdOk &= WriteReferenceEntry(refEntry);

        // then write any empty references remaining at end of file
        for ( int i = currentRefID+1; i < numReferences; ++i ) {
            BaiReferenceEntry emptyEntry(i);
            createdOk &= WriteReferenceEntry(emptyEntry);
        }
    }

    // rewind reader now that we're done building
    createdOk &= m_reader->Rewind();

    // return result
    return createdOk;
}

// returns format's file extension
const string BamStandardIndex::Extension(void) {
    return BamStandardIndex::BAI_EXTENSION;
}

bool BamStandardIndex::GetOffsets(const BamRegion& region, vector<int64_t>& offsets) {

    // cannot calculate offsets if unknown/invalid reference ID requested
    if ( region.LeftRefID < 0 || region.LeftRefID >= (int)m_indexFileSummary.size() )
        return false;

    // retrieve index summary for left bound reference
    const BaiReferenceSummary& refSummary = m_indexFileSummary.at(region.LeftRefID);

    // set up region boundaries based on actual BamReader data
    uint32_t begin;
    uint32_t end;
    if ( !AdjustRegion(region, begin, end) ) {
        cerr << "BamStandardIndex ERROR: cannot calculate offsets on invalid region" << endl;
        return false;
    }

    // retrieve all candidate bin IDs for region
    set<uint16_t> candidateBins;
    CalculateCandidateBins(begin, end, candidateBins);

    // use reference's linear offsets to calculate the minimum offset
    // that must be considered to find overlap
    const uint64_t& minOffset = CalculateMinOffset(refSummary, begin);

    // attempt to use reference summary, minOffset, & candidateBins to calculate offsets
    // no data should not be error
    if ( !CalculateCandidateOffsets(refSummary, minOffset, candidateBins, offsets) ) {
        cerr << "BamStandardIndex ERROR: could not calculate candidate offsets for requested region" << endl;
        return false;
    }

    // ensure that offsets are sorted before returning
    sort( offsets.begin(), offsets.end() );

    // return succes
    return true;
}

// returns whether reference has alignments or no
bool BamStandardIndex::HasAlignments(const int& referenceID) const {
    if ( referenceID < 0 || referenceID >= (int)m_indexFileSummary.size() )
        return false;
    const BaiReferenceSummary& refSummary = m_indexFileSummary.at(referenceID);
    return ( refSummary.NumBins > 0 );
}

bool BamStandardIndex::IsFileOpen(void) const {
    return ( m_indexStream != 0 );
}

// attempts to use index data to jump to @region, returns success/fail
// a "successful" jump indicates no error, but not whether this region has data
//   * thus, the method sets a flag to indicate whether there are alignments
//     available after the jump position
bool BamStandardIndex::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) {

    // clear out flag
    *hasAlignmentsInRegion = false;

    // skip if reader is not valid or is not open
    if ( m_reader == 0 || !m_reader->IsOpen() )
        return false;

    // calculate offsets for this region
    vector<int64_t> offsets;
    if ( !GetOffsets(region, offsets) ) {
        cerr << "BamStandardIndex ERROR: could not jump"
             << ", unable to retrieve offsets for region" << endl;
        return false;
    }

    // iterate through candidate offsets
    BamAlignment al;
    vector<int64_t>::const_iterator offsetIter = offsets.begin();
    vector<int64_t>::const_iterator offsetEnd  = offsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter) {

        // attempt seek
        if ( !m_reader->Seek(*offsetIter) ) {
            cerr << "BamStandardIndex ERROR: could not jump"
                 << ", there was a problem seeking in BAM file" << endl;
            return false;
        }

        // load first available alignment, setting flag to true if data exists
        *hasAlignmentsInRegion = m_reader->LoadNextAlignment(al);

        // if this alignment corresponds to desired position
        // return success of seeking back to the offset before the 'current offset' (to cover overlaps)
        if ( ((al.RefID == region.LeftRefID) &&
             ((al.Position + al.Length) > region.LeftPosition)) ||
             (al.RefID > region.LeftRefID) )
        {
            if ( offsetIter != offsets.begin() )
                --offsetIter;
            return m_reader->Seek(*offsetIter);
        }
    }

    // return success (no offset data is not an error,
    // but hasAlignments flag will be marked accordingly)
    return true;
}

// loads existing data from file into memory
bool BamStandardIndex::Load(const std::string& filename) {

    // attempt open index file (read-only)
    if ( !OpenFile(filename, "rb") ) {
        cerr << "BamStandardIndex ERROR: could not open input index file: " << filename
             << ", aborting index load" << endl;
        return false;
    }

    // if invalid format 'magic number', close & return failure
    if ( !CheckMagicNumber() ) {
        cerr << "BamStandardIndex ERROR: unexpected format for index file: " << filename
             << ", aborting index load" << endl;
        CloseFile();
        return false;
    }

    // attempt to load index file summary, return success/failure
    if ( !SummarizeIndexFile() ) {
        cerr << "BamStandardIndex ERROR: could not generate a summary of index file " << filename
             << ", aborting index load" << endl;
        CloseFile();
        return false;
    }

    // if we get here, index summary is loaded OK
    return true;
}

uint64_t BamStandardIndex::LookupLinearOffset(const BaiReferenceSummary& refSummary, const int& index) {

    // attempt seek to proper index file position
    const int64_t linearOffsetFilePosition = (int64_t)refSummary.FirstLinearOffsetFilePosition +
                                             index*BamStandardIndex::SIZEOF_LINEAROFFSET;
    if ( !Seek(linearOffsetFilePosition, SEEK_SET) )
        return 0;

    // read linear offset from BAI file
    uint64_t linearOffset(0);
    if ( !ReadLinearOffset(linearOffset) )
        return 0;
    return linearOffset;
}

void BamStandardIndex::MergeAlignmentChunks(BaiAlignmentChunkVector& chunks) {

    // skip if chunks are empty, nothing to merge
    if ( chunks.empty() )
        return;

    // set up merged alignment chunk container
    BaiAlignmentChunkVector mergedChunks;
    mergedChunks.push_back( chunks[0] );

    // iterate over chunks
    int i = 0;
    BaiAlignmentChunkVector::iterator chunkIter = chunks.begin();
    BaiAlignmentChunkVector::iterator chunkEnd  = chunks.end();
    for ( ++chunkIter; chunkIter != chunkEnd; ++chunkIter) {

        // get 'currentMergeChunk' based on numeric index
        BaiAlignmentChunk& currentMergeChunk = mergedChunks[i];

        // get sourceChunk based on source vector iterator
        BaiAlignmentChunk& sourceChunk = (*chunkIter);

        // if currentMergeChunk ends where sourceChunk starts, then merge the two
        if ( currentMergeChunk.Stop>>16 == sourceChunk.Start>>16 )
            currentMergeChunk.Stop = sourceChunk.Stop;

        // otherwise
        else {
            // append sourceChunk after currentMergeChunk
            mergedChunks.push_back(sourceChunk);

            // update i, so the next iteration will consider the
            // recently-appended sourceChunk as new mergeChunk candidate
            ++i;
        }
    }

    // saved newly-merged chunks into (parameter) chunks
    chunks = mergedChunks;
}

bool BamStandardIndex::OpenFile(const std::string& filename, const char* mode) {

    // make sure any previous index file is closed
    CloseFile();

    // attempt to open file
    m_indexStream = fopen(filename.c_str(), mode);
    return IsFileOpen();
}

bool BamStandardIndex::ReadBinID(uint32_t& binId) {
    size_t elementsRead = 0;
    elementsRead += fread(&binId, sizeof(binId), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(binId);
    return ( elementsRead == 1 );
}

bool BamStandardIndex::ReadBinIntoBuffer(uint32_t& binId, int32_t& numAlignmentChunks) {

    bool readOk = true;

    // read bin header
    readOk &= ReadBinID(binId);
    readOk &= ReadNumAlignmentChunks(numAlignmentChunks);

    // read bin contents
    const unsigned int bytesRequested = numAlignmentChunks*BamStandardIndex::SIZEOF_ALIGNMENTCHUNK;
    readOk &= ReadIntoBuffer(bytesRequested);

    // return success/failure
    return readOk;
}

bool BamStandardIndex::ReadIntoBuffer(const unsigned int& bytesRequested) {

    // ensure that our buffer is big enough for request
    BamStandardIndex::CheckBufferSize(m_buffer, m_bufferLength, bytesRequested);

    // read from BAI file stream
    size_t bytesRead = fread( m_buffer, sizeof(char), bytesRequested, m_indexStream );
    return ( bytesRead == (size_t)bytesRequested );
}

bool BamStandardIndex::ReadLinearOffset(uint64_t& linearOffset) {
    size_t elementsRead = 0;
    elementsRead += fread(&linearOffset, sizeof(linearOffset), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_64(linearOffset);
    return ( elementsRead == 1 );
}

bool BamStandardIndex::ReadNumAlignmentChunks(int& numAlignmentChunks) {
    size_t elementsRead = 0;
    elementsRead += fread(&numAlignmentChunks, sizeof(numAlignmentChunks), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numAlignmentChunks);
    return ( elementsRead == 1 );
}

bool BamStandardIndex::ReadNumBins(int& numBins) {
    size_t elementsRead = 0;
    elementsRead += fread(&numBins, sizeof(numBins), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numBins);
    return ( elementsRead == 1 );
}

bool BamStandardIndex::ReadNumLinearOffsets(int& numLinearOffsets) {
    size_t elementsRead = 0;
    elementsRead += fread(&numLinearOffsets, sizeof(numLinearOffsets), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numLinearOffsets);
    return ( elementsRead == 1 );
}

bool BamStandardIndex::ReadNumReferences(int& numReferences) {
    size_t elementsRead = 0;
    elementsRead += fread(&numReferences, sizeof(numReferences), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    return ( elementsRead == 1 );
}

void BamStandardIndex::ReserveForSummary(const int& numReferences) {
    m_indexFileSummary.clear();
    m_indexFileSummary.assign( numReferences, BaiReferenceSummary() );
}

void BamStandardIndex::SaveAlignmentChunkToBin(BaiBinMap& binMap,
                                    const uint32_t& currentBin,
                                    const uint64_t& currentOffset,
                                    const uint64_t& lastOffset)
{
    // create new alignment chunk
    BaiAlignmentChunk newChunk(currentOffset, lastOffset);



    // if no entry exists yet for this bin, create one and store alignment chunk
    BaiBinMap::iterator binIter = binMap.find(currentBin);
    if ( binIter == binMap.end() ) {
        BaiAlignmentChunkVector newChunks;
        newChunks.push_back(newChunk);
        binMap.insert( pair<uint32_t, BaiAlignmentChunkVector>(currentBin, newChunks));
    }

    // otherwise, just append alignment chunk
    else {
        BaiAlignmentChunkVector& binChunks = (*binIter).second;
        binChunks.push_back( newChunk );
    }
}

void BamStandardIndex::SaveBinsSummary(const int& refId, const int& numBins) {
    BaiReferenceSummary& refSummary = m_indexFileSummary.at(refId);
    refSummary.NumBins = numBins;
    refSummary.FirstBinFilePosition = Tell();
}

void BamStandardIndex::SaveLinearOffsetEntry(BaiLinearOffsetVector& offsets,
                                             const int& alignmentStartPosition,
                                             const int& alignmentStopPosition,
                                             const uint64_t& lastOffset)
{
    // get converted offsets
    const int beginOffset = alignmentStartPosition >> BamStandardIndex::BAM_LIDX_SHIFT;
    const int endOffset   = (alignmentStopPosition - 1) >> BamStandardIndex::BAM_LIDX_SHIFT;

    // resize vector if necessary
    int oldSize = offsets.size();
    int newSize = endOffset + 1;
    if ( oldSize < newSize )
        offsets.resize(newSize, 0);

    // store offset
    for( int i = beginOffset + 1; i <= endOffset; ++i ) {
        if ( offsets[i] == 0 )
            offsets[i] = lastOffset;
    }
}

void BamStandardIndex::SaveLinearOffsetsSummary(const int& refId, const int& numLinearOffsets) {
    BaiReferenceSummary& refSummary = m_indexFileSummary.at(refId);
    refSummary.NumLinearOffsets = numLinearOffsets;
    refSummary.FirstLinearOffsetFilePosition = Tell();
}

// seek to position in index file stream
bool BamStandardIndex::Seek(const int64_t& position, const int& origin) {
    return ( fseek64(m_indexStream, position, origin) == 0 );
}

// change the index caching behavior
void BamStandardIndex::SetCacheMode(const BamIndex::IndexCacheMode& mode) {
    m_cacheMode = mode;
    // do nothing else here ? cache mode will be ignored from now on, most likely
}

bool BamStandardIndex::SkipBins(const int& numBins) {
    uint32_t binId;
    int32_t numAlignmentChunks;
    bool skippedOk = true;
    for (int i = 0; i < numBins; ++i)
        skippedOk &= ReadBinIntoBuffer(binId, numAlignmentChunks); // results & buffer ignored
    return skippedOk;
}

bool BamStandardIndex::SkipLinearOffsets(const int& numLinearOffsets) {
    const unsigned int bytesRequested = numLinearOffsets*BamStandardIndex::SIZEOF_LINEAROFFSET;
    return ReadIntoBuffer(bytesRequested);
}

void BamStandardIndex::SortLinearOffsets(BaiLinearOffsetVector& linearOffsets) {
    sort( linearOffsets.begin(), linearOffsets.end() );
}

bool BamStandardIndex::SummarizeBins(BaiReferenceSummary& refSummary) {

    // load number of bins
    int numBins;
    if ( !ReadNumBins(numBins) )
        return false;

    // store bins summary for this reference
    refSummary.NumBins = numBins;
    refSummary.FirstBinFilePosition = Tell();

    // attempt skip reference bins, return success/failure
    if ( !SkipBins(numBins) )
        return false;

    // if we get here, bin summarized OK
    return true;
}

bool BamStandardIndex::SummarizeIndexFile(void) {

    // load number of reference sequences
    int numReferences;
    if ( !ReadNumReferences(numReferences) )
        return false;

    // initialize file summary data
    ReserveForSummary(numReferences);

    // iterate over reference entries
    bool loadedOk = true;
    BaiFileSummary::iterator summaryIter = m_indexFileSummary.begin();
    BaiFileSummary::iterator summaryEnd  = m_indexFileSummary.end();
    for ( int i = 0; summaryIter != summaryEnd; ++summaryIter, ++i )
        loadedOk &= SummarizeReference(*summaryIter);

    // return result
    return loadedOk;
}

bool BamStandardIndex::SummarizeLinearOffsets(BaiReferenceSummary& refSummary) {

    // load number of linear offsets
    int numLinearOffsets;
    if ( !ReadNumLinearOffsets(numLinearOffsets) )
        return false;

    // store bin summary data for this reference
    refSummary.NumLinearOffsets = numLinearOffsets;
    refSummary.FirstLinearOffsetFilePosition = Tell();

    // skip linear offsets in index file
    if ( !SkipLinearOffsets(numLinearOffsets) )
        return false;

    // if get here, linear offsets summarized OK
    return true;
}

bool BamStandardIndex::SummarizeReference(BaiReferenceSummary& refSummary) {

    bool loadedOk = true;
    loadedOk &= SummarizeBins(refSummary);
    loadedOk &= SummarizeLinearOffsets(refSummary);
    return loadedOk;
}

// return position of file pointer in index file stream
int64_t BamStandardIndex::Tell(void) const {
    return ftell64(m_indexStream);
}

bool BamStandardIndex::WriteAlignmentChunk(const BaiAlignmentChunk& chunk) {

    size_t elementsWritten = 0;

    // localize alignment chunk offsets
    uint64_t start = chunk.Start;
    uint64_t stop  = chunk.Stop;

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_64(start);
        SwapEndian_64(stop);
    }

    // write to index file
    elementsWritten += fwrite(&start, sizeof(start), 1, m_indexStream);
    elementsWritten += fwrite(&stop,  sizeof(stop),  1, m_indexStream);

    // return success/failure of write
    return ( elementsWritten == 2 );
}

bool BamStandardIndex::WriteAlignmentChunks(BaiAlignmentChunkVector& chunks) {

    // make sure chunks are merged (simplified) before writing & saving summary
    MergeAlignmentChunks(chunks);

    size_t elementsWritten = 0;

    // write chunks
    int32_t chunkCount = chunks.size();
    if ( m_isBigEndian ) SwapEndian_32(chunkCount);
    elementsWritten += fwrite(&chunkCount, sizeof(chunkCount), 1, m_indexStream);

    // iterate over chunks
    bool chunksOk = true;
    BaiAlignmentChunkVector::const_iterator chunkIter = chunks.begin();
    BaiAlignmentChunkVector::const_iterator chunkEnd  = chunks.end();
    for ( ; chunkIter != chunkEnd; ++chunkIter )
        chunksOk &= WriteAlignmentChunk( (*chunkIter) );

    // return success/failure of write
    return ( (elementsWritten == 1) && chunksOk );
}

bool BamStandardIndex::WriteBin(const uint32_t& binId, BaiAlignmentChunkVector& chunks) {

    size_t elementsWritten = 0;

    // write BAM bin ID
    uint32_t binKey = binId;
    if ( m_isBigEndian ) SwapEndian_32(binKey);
    elementsWritten += fwrite(&binKey, sizeof(binKey), 1, m_indexStream);

    // write bin's alignment chunks
    bool chunksOk = WriteAlignmentChunks(chunks);

    // return success/failure of write
    return ( (elementsWritten == 1) && chunksOk );
}

bool BamStandardIndex::WriteBins(const int& refId, BaiBinMap& bins) {

    size_t elementsWritten = 0;

    // write number of bins
    int32_t binCount = bins.size();
    if ( m_isBigEndian ) SwapEndian_32(binCount);
    elementsWritten += fwrite(&binCount, sizeof(binCount), 1, m_indexStream);

    // save summary for reference's bins
    SaveBinsSummary(refId, bins.size());

    // iterate over bins
    bool binsOk = true;
    BaiBinMap::iterator binIter = bins.begin();
    BaiBinMap::iterator binEnd  = bins.end();
    for ( ; binIter != binEnd; ++binIter )
        binsOk &= WriteBin( (*binIter).first, (*binIter).second );

    // return success/failure of write
    return ( (elementsWritten == 1) && binsOk );
}

bool BamStandardIndex::WriteHeader(void) {

    size_t elementsWritten = 0;

    // write magic number
    elementsWritten += fwrite(BamStandardIndex::BAI_MAGIC, sizeof(char), 4, m_indexStream);

    // write number of reference sequences
    int32_t numReferences = m_indexFileSummary.size();
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    elementsWritten += fwrite(&numReferences, sizeof(numReferences), 1, m_indexStream);

    // return success/failure of write
    return (elementsWritten == 5);
}

bool BamStandardIndex::WriteLinearOffsets(const int& refId, BaiLinearOffsetVector& linearOffsets) {

    // make sure linear offsets are sorted before writing & saving summary
    SortLinearOffsets(linearOffsets);

    size_t elementsWritten = 0;

    // write number of linear offsets
    int32_t offsetCount = linearOffsets.size();
    if ( m_isBigEndian ) SwapEndian_32(offsetCount);
    elementsWritten += fwrite(&offsetCount, sizeof(offsetCount), 1, m_indexStream);

    // save summary for reference's linear offsets
    SaveLinearOffsetsSummary(refId, linearOffsets.size());

    // iterate over linear offsets
    BaiLinearOffsetVector::const_iterator offsetIter = linearOffsets.begin();
    BaiLinearOffsetVector::const_iterator offsetEnd  = linearOffsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter ) {

        // write linear offset
        uint64_t linearOffset = (*offsetIter);
        if ( m_isBigEndian ) SwapEndian_64(linearOffset);
        elementsWritten += fwrite(&linearOffset, sizeof(linearOffset), 1, m_indexStream);
    }

    // return success/failure of write
    return ( elementsWritten == (size_t)(linearOffsets.size() + 1) );
}

bool BamStandardIndex::WriteReferenceEntry(BaiReferenceEntry& refEntry) {
    bool refOk = true;
    refOk &= WriteBins(refEntry.ID, refEntry.Bins);
    refOk &= WriteLinearOffsets(refEntry.ID, refEntry.LinearOffsets);
    return refOk;
}
