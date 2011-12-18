// ***************************************************************************
// BamStandardIndex.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the standardized BAM index format (".bai")
// ***************************************************************************

#include "api/BamAlignment.h"
#include "api/internal/bam/BamReader_p.h"
#include "api/internal/index/BamStandardIndex_p.h"
#include "api/internal/io/BamDeviceFactory_p.h"
#include "api/internal/utils/BamException_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <sstream>
using namespace std;

// -----------------------------------
// static BamStandardIndex constants
// -----------------------------------

const int BamStandardIndex::MAX_BIN               = 37450;  // =(8^6-1)/7+1
const int BamStandardIndex::BAM_LIDX_SHIFT        = 14;
const string BamStandardIndex::BAI_EXTENSION      = ".bai";
const char* const BamStandardIndex::BAI_MAGIC     = "BAI\1";
const int BamStandardIndex::SIZEOF_ALIGNMENTCHUNK = sizeof(uint64_t)*2;
const int BamStandardIndex::SIZEOF_BINCORE        = sizeof(uint32_t) + sizeof(int32_t);
const int BamStandardIndex::SIZEOF_LINEAROFFSET   = sizeof(uint64_t);

// ----------------------------
// RaiiWrapper implementation
// ----------------------------

BamStandardIndex::RaiiWrapper::RaiiWrapper(void)
    : Device(0)
    , Buffer(0)
{ }

BamStandardIndex::RaiiWrapper::~RaiiWrapper(void) {

    if ( Device ) {
        Device->Close();
        delete Device;
        Device = 0;
    }

    if ( Buffer ) {
        delete[] Buffer;
        Buffer = 0;
    }
}

// ---------------------------------
// BamStandardIndex implementation
// ---------------------------------

// ctor
BamStandardIndex::BamStandardIndex(Internal::BamReaderPrivate* reader)
    : BamIndex(reader)
    , m_bufferLength(0)
{
     m_isBigEndian = BamTools::SystemIsBigEndian();
}

// dtor
BamStandardIndex::~BamStandardIndex(void) {
    CloseFile();
}

void BamStandardIndex::AdjustRegion(const BamRegion& region, uint32_t& begin, uint32_t& end) {

    // retrieve references from reader
    const RefVector& references = m_reader->GetReferenceData();

    // LeftPosition cannot be greater than or equal to reference length
    if ( region.LeftPosition >= references.at(region.LeftRefID).RefLength )
        throw BamException("BamStandardIndex::AdjustRegion", "invalid region requested");

    // set region 'begin'
    begin = (unsigned int)region.LeftPosition;

    // if right bound specified AND left&right bounds are on same reference
    // OK to use right bound position as region 'end'
    if ( region.isRightBoundSpecified() && ( region.LeftRefID == region.RightRefID ) )
        end = (unsigned int)region.RightPosition;

    // otherwise, set region 'end' to last reference base
    else end = (unsigned int)references.at(region.LeftRefID).RefLength;
}

// [begin, end)
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

void BamStandardIndex::CalculateCandidateOffsets(const BaiReferenceSummary& refSummary,
                                                 const uint64_t& minOffset,
                                                 set<uint16_t>& candidateBins,
                                                 vector<int64_t>& offsets)
{
    // seek to first bin
    Seek(refSummary.FirstBinFilePosition, SEEK_SET);

    // iterate over reference bins
    uint32_t binId;
    int32_t numAlignmentChunks;
    set<uint16_t>::iterator candidateBinIter;
    for ( int i = 0; i < refSummary.NumBins; ++i ) {

        // read bin contents (if successful, alignment chunks are now in m_buffer)
        ReadBinIntoBuffer(binId, numAlignmentChunks);

        // see if bin is a 'candidate bin'
        candidateBinIter = candidateBins.find(binId);

        // if not, move on to next bin
        if ( candidateBinIter == candidateBins.end() )
            continue;

        // otherwise, check bin's contents against for overlap
        else {

            size_t offset = 0;
            uint64_t chunkStart;
            uint64_t chunkStop;

            // iterate over alignment chunks
            for ( int j = 0; j < numAlignmentChunks; ++j ) {

                // read chunk start & stop from buffer
                memcpy((char*)&chunkStart, m_resources.Buffer+offset, sizeof(uint64_t));
                offset += sizeof(uint64_t);
                memcpy((char*)&chunkStop, m_resources.Buffer+offset, sizeof(uint64_t));
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
    } catch ( std::bad_alloc&  ) {
        stringstream s("");
        s << "out of memory when allocating " << requestedBytes << " bytes";
        throw BamException("BamStandardIndex::CheckBufferSize", s.str());
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
    } catch ( std::bad_alloc& ) {
        stringstream s("");
        s << "out of memory when allocating " << requestedBytes << " bytes";
        throw BamException("BamStandardIndex::CheckBufferSize", s.str());
    }
}

void BamStandardIndex::CheckMagicNumber(void) {

    // check 'magic number' to see if file is BAI index
    char magic[4];
    const int64_t numBytesRead = m_resources.Device->Read(magic, sizeof(magic));
    if ( numBytesRead != 4 )
        throw BamException("BamStandardIndex::CheckMagicNumber", "could not read BAI magic number");

    // compare to expected value
    if ( strncmp(magic, BamStandardIndex::BAI_MAGIC, 4) != 0 )
        throw BamException("BamStandardIndex::CheckMagicNumber", "invalid BAI magic number");
}

void BamStandardIndex::ClearReferenceEntry(BaiReferenceEntry& refEntry) {
    refEntry.ID = -1;
    refEntry.Bins.clear();
    refEntry.LinearOffsets.clear();
}

void BamStandardIndex::CloseFile(void) {

    // close file stream
    if ( IsDeviceOpen() ) {
        m_resources.Device->Close();
        delete m_resources.Device;
        m_resources.Device = 0;
    }

    // clear index file summary data
    m_indexFileSummary.clear();

    // clean up I/O buffer
    delete[] m_resources.Buffer;
    m_resources.Buffer = 0;
    m_bufferLength = 0;
}

// builds index from associated BAM file & writes out to index file
bool BamStandardIndex::Create(void) {

    // skip if BamReader is invalid or not open
    if ( m_reader == 0 || !m_reader->IsOpen() ) {
        SetErrorString("BamStandardIndex::Create", "could not create index: reader is not open");
        return false;
    }

    // rewind BamReader
    if ( !m_reader->Rewind() ) {
        const string readerError = m_reader->GetErrorString();
        const string message = "could not create index: \n\t" + readerError;
        SetErrorString("BamStandardIndex::Create", message);
        return false;
    }

    try {

        // open new index file (read & write)
        string indexFilename = m_reader->Filename() + Extension();
        OpenFile(indexFilename, IBamIODevice::ReadWrite);

        // initialize BaiFileSummary with number of references
        const int& numReferences = m_reader->GetReferenceCount();
        ReserveForSummary(numReferences);

        // initialize output file
        WriteHeader();

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
                    WriteReferenceEntry(refEntry);
                    ClearReferenceEntry(refEntry);

                    // write any empty references between (but *NOT* including) lastRefID & al.RefID
                    for ( int i = lastRefID+1; i < al.RefID; ++i ) {
                        BaiReferenceEntry emptyEntry(i);
                        WriteReferenceEntry(emptyEntry);
                    }

                    // update bin markers
                    currentOffset = lastOffset;
                    currentBin    = al.Bin;
                    lastBin       = al.Bin;
                    currentRefID  = al.RefID;
                }

                // otherwise, this is first pass
                // be sure to write any empty references up to (but *NOT* including) current RefID
                else {
                    for ( int i = 0; i < al.RefID; ++i ) {
                        BaiReferenceEntry emptyEntry(i);
                        WriteReferenceEntry(emptyEntry);
                    }
                }

                // update reference markers
                refEntry.ID = al.RefID;
                lastRefID   = al.RefID;
                lastBin     = defaultValue;
            }

            // if lastPosition greater than current alignment position - file not sorted properly
            else if ( lastPosition > al.Position ) {
                stringstream s("");
                s << "BAM file is not properly sorted by coordinate" << endl
                  << "Current alignment position: " << al.Position
                  << " < previous alignment position: " << lastPosition
                  << " on reference ID: " << al.RefID << endl;
                SetErrorString("BamStandardIndex::Create", s.str());
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
                SetErrorString("BamStandardIndex::Create", "calculating offsets failed");
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
            WriteReferenceEntry(refEntry);

            // then write any empty references remaining at end of file
            for ( int i = currentRefID+1; i < numReferences; ++i ) {
                BaiReferenceEntry emptyEntry(i);
                WriteReferenceEntry(emptyEntry);
            }
        }

    } catch ( BamException& e) {
        m_errorString = e.what();
        return false;
    }

    // rewind BamReader
    if ( !m_reader->Rewind() ) {
        const string readerError = m_reader->GetErrorString();
        const string message = "could not create index: \n\t" + readerError;
        SetErrorString("BamStandardIndex::Create", message);
        return false;
    }

    // return success
    return true;
}

// returns format's file extension
const string BamStandardIndex::Extension(void) {
    return BamStandardIndex::BAI_EXTENSION;
}

void BamStandardIndex::GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion) {

    // cannot calculate offsets if unknown/invalid reference ID requested
    if ( region.LeftRefID < 0 || region.LeftRefID >= (int)m_indexFileSummary.size() )
        throw BamException("BamStandardIndex::GetOffset", "invalid reference ID requested");

    // retrieve index summary for left bound reference
    const BaiReferenceSummary& refSummary = m_indexFileSummary.at(region.LeftRefID);

    // set up region boundaries based on actual BamReader data
    uint32_t begin;
    uint32_t end;
    AdjustRegion(region, begin, end);

    // retrieve all candidate bin IDs for region
    set<uint16_t> candidateBins;
    CalculateCandidateBins(begin, end, candidateBins);

    // use reference's linear offsets to calculate the minimum offset
    // that must be considered to find overlap
    const uint64_t& minOffset = CalculateMinOffset(refSummary, begin);

    // attempt to use reference summary, minOffset, & candidateBins to calculate offsets
    // no data should not be error, just bail
    vector<int64_t> offsets;
    CalculateCandidateOffsets(refSummary, minOffset, candidateBins, offsets);
    if ( offsets.empty() )
        return;
    
    // ensure that offsets are sorted before processing
    sort( offsets.begin(), offsets.end() );

    // binary search for an overlapping block (may not be first one though)
    BamAlignment al;
    typedef vector<int64_t>::const_iterator OffsetConstIterator;
    OffsetConstIterator offsetFirst = offsets.begin();
    OffsetConstIterator offsetIter  = offsetFirst;
    OffsetConstIterator offsetLast  = offsets.end();
    iterator_traits<OffsetConstIterator>::difference_type count = distance(offsetFirst, offsetLast);
    iterator_traits<OffsetConstIterator>::difference_type step;
    while ( count > 0 ) {
        offsetIter = offsetFirst;
        step = count/2;
        advance(offsetIter, step);

        // attempt seek to candidate offset
        const int64_t& candidateOffset = (*offsetIter);
        if ( !m_reader->Seek(candidateOffset) ) {
            const string readerError = m_reader->GetErrorString();
            const string message = "could not seek in BAM file: \n\t" + readerError;
            throw BamException("BamToolsIndex::GetOffset", message);
        }

        // load first available alignment, setting flag to true if data exists
        *hasAlignmentsInRegion = m_reader->LoadNextAlignment(al);

        // check alignment against region
        if ( al.GetEndPosition() <= region.LeftPosition ) {
            offsetFirst = ++offsetIter;
            count -= step+1;
        } else count = step;
    }

    // step back to the offset before the 'current offset' (to make sure we cover overlaps)
    if ( offsetIter != offsets.begin() )
        --offsetIter;
    offset = (*offsetIter);
}

// returns whether reference has alignments or no
bool BamStandardIndex::HasAlignments(const int& referenceID) const {
    if ( referenceID < 0 || referenceID >= (int)m_indexFileSummary.size() )
        return false;
    const BaiReferenceSummary& refSummary = m_indexFileSummary.at(referenceID);
    return ( refSummary.NumBins > 0 );
}

bool BamStandardIndex::IsDeviceOpen(void) const {
    if ( m_resources.Device == 0 )
        return false;
    return m_resources.Device->IsOpen();
}

// attempts to use index data to jump to @region, returns success/fail
// a "successful" jump indicates no error, but not whether this region has data
//   * thus, the method sets a flag to indicate whether there are alignments
//     available after the jump position
bool BamStandardIndex::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) {

    // clear out flag
    *hasAlignmentsInRegion = false;

    // skip if invalid reader or not open
    if ( m_reader == 0 || !m_reader->IsOpen() ) {
        SetErrorString("BamStandardIndex::Jump", "could not jump: reader is not open");
        return false;
    }

    // calculate nearest offset to jump to
    int64_t offset;
    try {
        GetOffset(region, offset, hasAlignmentsInRegion);
    } catch ( BamException& e ) {
        m_errorString = e.what();
        return false;
    }

    // if region has alignments, return success/fail of seeking there
    if ( *hasAlignmentsInRegion )
        return m_reader->Seek(offset);

    // otherwise, simply return true (but hasAlignmentsInRegion flag has been set to false)
    // (this is OK, BamReader will check this flag before trying to load data)
    return true;
}

// loads existing data from file into memory
bool BamStandardIndex::Load(const std::string& filename) {

    try {

        // attempt to open file (read-only)
        OpenFile(filename, IBamIODevice::ReadOnly);

        // validate format
        CheckMagicNumber();

        // load in-memory summary of index data
        SummarizeIndexFile();

        // return success
        return true;

    } catch ( BamException& e ) {
        m_errorString = e.what();
        return false;
    }
}

uint64_t BamStandardIndex::LookupLinearOffset(const BaiReferenceSummary& refSummary, const int& index) {

    // attempt seek to proper index file position
    const int64_t linearOffsetFilePosition = (int64_t)refSummary.FirstLinearOffsetFilePosition +
                                             index*BamStandardIndex::SIZEOF_LINEAROFFSET;
    Seek(linearOffsetFilePosition, SEEK_SET);

    // read linear offset from BAI file
    uint64_t linearOffset;
    ReadLinearOffset(linearOffset);
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

void BamStandardIndex::OpenFile(const std::string& filename, IBamIODevice::OpenMode mode) {

    // make sure any previous index file is closed
    CloseFile();

    m_resources.Device = BamDeviceFactory::CreateDevice(filename);
    if ( m_resources.Device == 0 ) {
        const string message = string("could not open file: ") + filename;
        throw BamException("BamStandardIndex::OpenFile", message);
    }

    // attempt to open file
    m_resources.Device->Open(mode);
    if ( !IsDeviceOpen() ) {
        const string message = string("could not open file: ") + filename;
        throw BamException("BamStandardIndex::OpenFile", message);
    }
}

void BamStandardIndex::ReadBinID(uint32_t& binId) {
    const int64_t numBytesRead = m_resources.Device->Read((char*)&binId, sizeof(binId));
    if ( m_isBigEndian ) SwapEndian_32(binId);
    if ( numBytesRead != sizeof(binId) )
        throw BamException("BamStandardIndex::ReadBinID", "could not read BAI bin ID");
}

void BamStandardIndex::ReadBinIntoBuffer(uint32_t& binId, int32_t& numAlignmentChunks) {

    // read bin header
    ReadBinID(binId);
    ReadNumAlignmentChunks(numAlignmentChunks);

    // read bin contents
    const unsigned int bytesRequested = numAlignmentChunks*BamStandardIndex::SIZEOF_ALIGNMENTCHUNK;
    ReadIntoBuffer(bytesRequested);
}

void BamStandardIndex::ReadIntoBuffer(const unsigned int& bytesRequested) {

    // ensure that our buffer is big enough for request
    BamStandardIndex::CheckBufferSize(m_resources.Buffer, m_bufferLength, bytesRequested);

    // read from BAI file stream
    const int64_t bytesRead = m_resources.Device->Read(m_resources.Buffer, bytesRequested);
    if ( bytesRead != (int64_t)bytesRequested ) {
        stringstream s("");
        s << "expected to read: " << bytesRequested << " bytes, "
          << "but instead read: " << bytesRead;
        throw BamException("BamStandardIndex::ReadIntoBuffer", s.str());
    }
}

void BamStandardIndex::ReadLinearOffset(uint64_t& linearOffset) {
    const int64_t numBytesRead = m_resources.Device->Read((char*)&linearOffset, sizeof(linearOffset));
    if ( m_isBigEndian ) SwapEndian_64(linearOffset);
    if ( numBytesRead != sizeof(linearOffset) )
        throw BamException("BamStandardIndex::ReadLinearOffset", "could not read BAI linear offset");
}

void BamStandardIndex::ReadNumAlignmentChunks(int& numAlignmentChunks) {
    const int64_t numBytesRead = m_resources.Device->Read((char*)&numAlignmentChunks, sizeof(numAlignmentChunks));
    if ( m_isBigEndian ) SwapEndian_32(numAlignmentChunks);
    if ( numBytesRead != sizeof(numAlignmentChunks) )
        throw BamException("BamStandardIndex::ReadNumAlignmentChunks", "could not read BAI chunk count");
}

void BamStandardIndex::ReadNumBins(int& numBins) {
    const int64_t numBytesRead = m_resources.Device->Read((char*)&numBins, sizeof(numBins));
    if ( m_isBigEndian ) SwapEndian_32(numBins);
    if ( numBytesRead != sizeof(numBins) )
        throw BamException("BamStandardIndex::ReadNumBins", "could not read BAI bin count");
}

void BamStandardIndex::ReadNumLinearOffsets(int& numLinearOffsets) {
    const int64_t numBytesRead = m_resources.Device->Read((char*)&numLinearOffsets, sizeof(numLinearOffsets));
    if ( m_isBigEndian ) SwapEndian_32(numLinearOffsets);
    if ( numBytesRead != sizeof(numLinearOffsets) )
        throw BamException("BamStandardIndex::ReadNumAlignmentChunks", "could not read BAI linear offset count");
}

void BamStandardIndex::ReadNumReferences(int& numReferences) {
    const int64_t numBytesRead = m_resources.Device->Read((char*)&numReferences, sizeof(numReferences));
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    if ( numBytesRead != sizeof(numReferences) )
        throw BamException("BamStandardIndex::ReadNumReferences", "could not read reference count");
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
void BamStandardIndex::Seek(const int64_t& position, const int origin) {
    if ( !m_resources.Device->Seek(position, origin) )
        throw BamException("BamStandardIndex::Seek", "could not seek in BAI file");
}

void BamStandardIndex::SkipBins(const int& numBins) {
    uint32_t binId;
    int32_t numAlignmentChunks;
    for (int i = 0; i < numBins; ++i)
        ReadBinIntoBuffer(binId, numAlignmentChunks); // results & buffer ignored
}

void BamStandardIndex::SkipLinearOffsets(const int& numLinearOffsets) {
    const unsigned int bytesRequested = numLinearOffsets*BamStandardIndex::SIZEOF_LINEAROFFSET;
    ReadIntoBuffer(bytesRequested);
}

void BamStandardIndex::SortLinearOffsets(BaiLinearOffsetVector& linearOffsets) {
    sort( linearOffsets.begin(), linearOffsets.end() );
}

void BamStandardIndex::SummarizeBins(BaiReferenceSummary& refSummary) {

    // load number of bins
    int numBins;
    ReadNumBins(numBins);

    // store bins summary for this reference
    refSummary.NumBins = numBins;
    refSummary.FirstBinFilePosition = Tell();

    // skip this reference's bins
    SkipBins(numBins);
}

void BamStandardIndex::SummarizeIndexFile(void) {

    // load number of reference sequences
    int numReferences;
    ReadNumReferences(numReferences);

    // initialize file summary data
    ReserveForSummary(numReferences);

    // iterate over reference entries
    BaiFileSummary::iterator summaryIter = m_indexFileSummary.begin();
    BaiFileSummary::iterator summaryEnd  = m_indexFileSummary.end();
    for ( int i = 0; summaryIter != summaryEnd; ++summaryIter, ++i )
        SummarizeReference(*summaryIter);
}

void BamStandardIndex::SummarizeLinearOffsets(BaiReferenceSummary& refSummary) {

    // load number of linear offsets
    int numLinearOffsets;
    ReadNumLinearOffsets(numLinearOffsets);

    // store bin summary data for this reference
    refSummary.NumLinearOffsets = numLinearOffsets;
    refSummary.FirstLinearOffsetFilePosition = Tell();

    // skip linear offsets in index file
    SkipLinearOffsets(numLinearOffsets);
}

void BamStandardIndex::SummarizeReference(BaiReferenceSummary& refSummary) {
    SummarizeBins(refSummary);
    SummarizeLinearOffsets(refSummary);
}

// return position of file pointer in index file stream
int64_t BamStandardIndex::Tell(void) const {
    return m_resources.Device->Tell();
}

void BamStandardIndex::WriteAlignmentChunk(const BaiAlignmentChunk& chunk) {

    // localize alignment chunk offsets
    uint64_t start = chunk.Start;
    uint64_t stop  = chunk.Stop;

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_64(start);
        SwapEndian_64(stop);
    }

    // write to index file
    int64_t numBytesWritten = 0;
    numBytesWritten += m_resources.Device->Write((const char*)&start, sizeof(start));
    numBytesWritten += m_resources.Device->Write((const char*)&stop, sizeof(stop));
    if ( numBytesWritten != (sizeof(start)+sizeof(stop)) )
        throw BamException("BamStandardIndex::WriteAlignmentChunk", "could not write BAI alignment chunk");
}

void BamStandardIndex::WriteAlignmentChunks(BaiAlignmentChunkVector& chunks) {

    // make sure chunks are merged (simplified) before writing & saving summary
    MergeAlignmentChunks(chunks);

    // write chunks
    int32_t chunkCount = chunks.size();
    if ( m_isBigEndian ) SwapEndian_32(chunkCount);
    const int64_t numBytesWritten = m_resources.Device->Write((const char*)&chunkCount, sizeof(chunkCount));
    if ( numBytesWritten != sizeof(chunkCount) )
        throw BamException("BamStandardIndex::WriteAlignmentChunks", "could not write BAI chunk count");

    // iterate over chunks
    BaiAlignmentChunkVector::const_iterator chunkIter = chunks.begin();
    BaiAlignmentChunkVector::const_iterator chunkEnd  = chunks.end();
    for ( ; chunkIter != chunkEnd; ++chunkIter )
        WriteAlignmentChunk( (*chunkIter) );
}

void BamStandardIndex::WriteBin(const uint32_t& binId, BaiAlignmentChunkVector& chunks) {

    // write BAM bin ID
    uint32_t binKey = binId;
    if ( m_isBigEndian ) SwapEndian_32(binKey);
    const int64_t numBytesWritten = m_resources.Device->Write((const char*)&binKey, sizeof(binKey));
    if ( numBytesWritten != sizeof(binKey) )
        throw BamException("BamStandardIndex::WriteBin", "could not write bin ID");

    // write bin's alignment chunks
    WriteAlignmentChunks(chunks);
}

void BamStandardIndex::WriteBins(const int& refId, BaiBinMap& bins) {

    // write number of bins
    int32_t binCount = bins.size();
    if ( m_isBigEndian ) SwapEndian_32(binCount);
    const int64_t numBytesWritten = m_resources.Device->Write((const char*)&binCount, sizeof(binCount));
    if ( numBytesWritten != sizeof(binCount) )
        throw BamException("BamStandardIndex::WriteBins", "could not write bin count");

    // save summary for reference's bins
    SaveBinsSummary(refId, bins.size());

    // iterate over bins
    BaiBinMap::iterator binIter = bins.begin();
    BaiBinMap::iterator binEnd  = bins.end();
    for ( ; binIter != binEnd; ++binIter )
        WriteBin( (*binIter).first, (*binIter).second );
}

void BamStandardIndex::WriteHeader(void) {

    int64_t numBytesWritten = 0;

    // write magic number
    numBytesWritten += m_resources.Device->Write(BamStandardIndex::BAI_MAGIC, 4);

    // write number of reference sequences
    int32_t numReferences = m_indexFileSummary.size();
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    numBytesWritten += m_resources.Device->Write((const char*)&numReferences, sizeof(numReferences));

    if ( numBytesWritten != sizeof(numReferences)+4 )
        throw BamException("BamStandardIndex::WriteHeader", "could not write BAI header");
}

void BamStandardIndex::WriteLinearOffsets(const int& refId, BaiLinearOffsetVector& linearOffsets) {

    // make sure linear offsets are sorted before writing & saving summary
    SortLinearOffsets(linearOffsets);

    int64_t numBytesWritten = 0;

    // write number of linear offsets
    int32_t offsetCount = linearOffsets.size();
    if ( m_isBigEndian ) SwapEndian_32(offsetCount);
    numBytesWritten += m_resources.Device->Write((const char*)&offsetCount, sizeof(offsetCount));

    // save summary for reference's linear offsets
    SaveLinearOffsetsSummary(refId, linearOffsets.size());

    // iterate over linear offsets
    BaiLinearOffsetVector::const_iterator offsetIter = linearOffsets.begin();
    BaiLinearOffsetVector::const_iterator offsetEnd  = linearOffsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter ) {

        // write linear offset
        uint64_t linearOffset = (*offsetIter);
        if ( m_isBigEndian ) SwapEndian_64(linearOffset);
        numBytesWritten += m_resources.Device->Write((const char*)&linearOffset, sizeof(linearOffset));
    }

    if ( numBytesWritten != (sizeof(offsetCount) + linearOffsets.size()*sizeof(uint64_t)) )
        throw BamException("BamStandardIndex::WriteLinearOffsets", "could not write BAI linear offsets");
}

void BamStandardIndex::WriteReferenceEntry(BaiReferenceEntry& refEntry) {
    WriteBins(refEntry.ID, refEntry.Bins);
    WriteLinearOffsets(refEntry.ID, refEntry.LinearOffsets);
}
