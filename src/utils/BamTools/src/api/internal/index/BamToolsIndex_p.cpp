// ***************************************************************************
// BamToolsIndex.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the BamTools index format (".bti")
// ***************************************************************************

#include "api/BamAlignment.h"
#include "api/internal/bam/BamReader_p.h"
#include "api/internal/index/BamToolsIndex_p.h"
#include "api/internal/io/BamDeviceFactory_p.h"
#include "api/internal/io/BgzfStream_p.h"
#include "api/internal/utils/BamException_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
using namespace std;

// --------------------------------
// static BamToolsIndex constants
// --------------------------------

const uint32_t BamToolsIndex::DEFAULT_BLOCK_LENGTH = 1000;
const string BamToolsIndex::BTI_EXTENSION     = ".bti";
const char* const BamToolsIndex::BTI_MAGIC    = "BTI\1";
const int BamToolsIndex::SIZEOF_BLOCK         = sizeof(int32_t)*2 + sizeof(int64_t);

// ----------------------------
// RaiiWrapper implementation
// ----------------------------

BamToolsIndex::RaiiWrapper::RaiiWrapper(void)
    : Device(0)
{ }

BamToolsIndex::RaiiWrapper::~RaiiWrapper(void) {
    if ( Device ) {
        Device->Close();
        delete Device;
        Device = 0;
    }
}

// ------------------------------
// BamToolsIndex implementation
// ------------------------------

// ctor
BamToolsIndex::BamToolsIndex(Internal::BamReaderPrivate* reader)
    : BamIndex(reader)
    , m_blockSize(BamToolsIndex::DEFAULT_BLOCK_LENGTH)
    , m_inputVersion(0)
    , m_outputVersion(BTI_2_0) // latest version - used for writing new index files
{
    m_isBigEndian = BamTools::SystemIsBigEndian();
}

// dtor
BamToolsIndex::~BamToolsIndex(void) {
    CloseFile();
}

void BamToolsIndex::CheckMagicNumber(void) {

    // read magic number
    char magic[4];
    const int64_t numBytesRead = m_resources.Device->Read(magic, 4);
    if ( numBytesRead != 4 )
        throw BamException("BamToolsIndex::CheckMagicNumber", "could not read BTI magic number");

    // validate expected magic number
    if ( strncmp(magic, BamToolsIndex::BTI_MAGIC, 4) != 0 )
        throw BamException("BamToolsIndex::CheckMagicNumber", "invalid BTI magic number");
}

// check index file version, return true if OK
void BamToolsIndex::CheckVersion(void) {

    // read version from file
    const int64_t numBytesRead = m_resources.Device->Read((char*)&m_inputVersion, sizeof(m_inputVersion));
    if ( numBytesRead != sizeof(m_inputVersion) )
        throw BamException("BamToolsIndex::CheckVersion", "could not read format version");
    if ( m_isBigEndian ) SwapEndian_32(m_inputVersion);

    // if version is negative, or zero
    if ( m_inputVersion <= 0 )
        throw BamException("BamToolsIndex::CheckVersion", "invalid format version");

    // if version is newer than can be supported by this version of bamtools
    else if ( m_inputVersion > m_outputVersion ) {
        const string message = "unsupported format: this index was created by a newer version of BamTools. "
                               "Update your local version of BamTools to use the index file.";
        throw BamException("BamToolsIndex::CheckVersion", message);
    }

    // ------------------------------------------------------------------
    // check for deprecated, unsupported versions
    // (the format had to be modified to accomodate a particular bug fix)

    // Version 2.0: introduced support for half-open intervals, instead of the old closed intervals
    //   respondBy: throwing exception - we're not going to try to handle the old BTI files.
    else if ( (Version)m_inputVersion < BamToolsIndex::BTI_2_0 ) {
        const string message = "unsupported format: this version of the index may not properly handle "
                               "coordinate intervals. Please run 'bamtools index -bti -in yourData.bam' "
                               "to generate an up-to-date, fixed BTI file.";
        throw BamException("BamToolsIndex::CheckVersion", message);
    }
}

void BamToolsIndex::ClearReferenceEntry(BtiReferenceEntry& refEntry) {
    refEntry.ID = -1;
    refEntry.Blocks.clear();
}

void BamToolsIndex::CloseFile(void) {
    if ( IsDeviceOpen() ) {
        m_resources.Device->Close();
        delete m_resources.Device;
        m_resources.Device = 0;
    }
    m_indexFileSummary.clear();
}

// builds index from associated BAM file & writes out to index file
bool BamToolsIndex::Create(void) {

    // skip if BamReader is invalid or not open
    if ( m_reader == 0 || !m_reader->IsOpen() ) {
        SetErrorString("BamToolsIndex::Create", "could not create index: reader is not open");
        return false;
    }

    // rewind BamReader
    if ( !m_reader->Rewind() ) {
        const string readerError = m_reader->GetErrorString();
        const string message = "could not create index: \n\t" + readerError;
        SetErrorString("BamToolsIndex::Create", message);
        return false;
    }

    try {
        // open new index file (read & write)
        const string indexFilename = m_reader->Filename() + Extension();
        OpenFile(indexFilename, IBamIODevice::ReadWrite);

        // initialize BtiFileSummary with number of references
        const int& numReferences = m_reader->GetReferenceCount();
        InitializeFileSummary(numReferences);

        // intialize output file header
        WriteHeader();

        // index building markers
        uint32_t currentBlockCount      = 0;
        int64_t currentAlignmentOffset  = m_reader->Tell();
        int32_t blockRefId              = -1;
        int32_t blockMaxEndPosition     = -1;
        int64_t blockStartOffset        = currentAlignmentOffset;
        int32_t blockStartPosition      = -1;

        // plow through alignments, storing index entries
        BamAlignment al;
        BtiReferenceEntry refEntry;
        while ( m_reader->LoadNextAlignment(al) ) {

            // if moved to new reference
            if ( al.RefID != blockRefId ) {

                // if first pass, check:
                if ( currentBlockCount == 0 ) {

                    // write any empty references up to (but not including) al.RefID
                    for ( int i = 0; i < al.RefID; ++i )
                        WriteReferenceEntry( BtiReferenceEntry(i) );
                }

                // not first pass:
                else {

                    // store previous BTI block data in reference entry
                    const BtiBlock block(blockMaxEndPosition, blockStartOffset, blockStartPosition);
                    refEntry.Blocks.push_back(block);

                    // write reference entry, then clear
                    WriteReferenceEntry(refEntry);
                    ClearReferenceEntry(refEntry);

                    // write any empty references between (but not including)
                    // the last blockRefID and current al.RefID
                    for ( int i = blockRefId+1; i < al.RefID; ++i )
                        WriteReferenceEntry( BtiReferenceEntry(i) );

                    // reset block count
                    currentBlockCount = 0;
                }

                // set ID for new reference entry
                refEntry.ID = al.RefID;
            }

            // if beginning of block, update counters
            if ( currentBlockCount == 0 ) {
                blockRefId          = al.RefID;
                blockStartOffset    = currentAlignmentOffset;
                blockStartPosition  = al.Position;
                blockMaxEndPosition = al.GetEndPosition();
            }

            // increment block counter
            ++currentBlockCount;

            // check end position
            const int32_t alignmentEndPosition = al.GetEndPosition();
            if ( alignmentEndPosition > blockMaxEndPosition )
                blockMaxEndPosition = alignmentEndPosition;

            // if block is full, get offset for next block, reset currentBlockCount
            if ( currentBlockCount == m_blockSize ) {

                // store previous block data in reference entry
                const BtiBlock block(blockMaxEndPosition, blockStartOffset, blockStartPosition);
                refEntry.Blocks.push_back(block);

                // update markers
                blockStartOffset  = m_reader->Tell();
                currentBlockCount = 0;
            }

            // not the best name, but for the next iteration, this value will be the offset of the
            // *current* alignment. this is necessary because we won't know if this next alignment
            // is on a new reference until we actually read it
            currentAlignmentOffset = m_reader->Tell();
        }

        // after finishing alignments, if any data was read, check:
        if ( blockRefId >= 0 ) {

            // store last BTI block data in reference entry
            const BtiBlock block(blockMaxEndPosition, blockStartOffset, blockStartPosition);
            refEntry.Blocks.push_back(block);

            // write last reference entry, then clear
            WriteReferenceEntry(refEntry);
            ClearReferenceEntry(refEntry);

            // then write any empty references remaining at end of file
            for ( int i = blockRefId+1; i < numReferences; ++i )
                WriteReferenceEntry( BtiReferenceEntry(i) );
        }

    } catch ( BamException& e ) {
        m_errorString = e.what();
        return false;
    }

    // rewind BamReader
    if ( !m_reader->Rewind() ) {
        const string readerError = m_reader->GetErrorString();
        const string message = "could not create index: \n\t" + readerError;
        SetErrorString("BamToolsIndex::Create", message);
        return false;
    }

    // return success
    return true;
}

// returns format's file extension
const std::string BamToolsIndex::Extension(void) {
    return BamToolsIndex::BTI_EXTENSION;
}

void BamToolsIndex::GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion) {

    // return false ref ID is not a valid index in file summary data
    if ( region.LeftRefID < 0 || region.LeftRefID >= (int)m_indexFileSummary.size() )
        throw BamException("BamToolsIndex::GetOffset", "invalid region requested");

    // retrieve reference index data for left bound reference
    BtiReferenceEntry refEntry(region.LeftRefID);
    ReadReferenceEntry(refEntry);

    // binary search for an overlapping block (may not be first one though)
    bool found = false;
    typedef BtiBlockVector::const_iterator BtiBlockConstIterator;
    BtiBlockConstIterator blockFirst = refEntry.Blocks.begin();
    BtiBlockConstIterator blockIter  = blockFirst;
    BtiBlockConstIterator blockLast  = refEntry.Blocks.end();
    iterator_traits<BtiBlockConstIterator>::difference_type count = distance(blockFirst, blockLast);
    iterator_traits<BtiBlockConstIterator>::difference_type step;
    while ( count > 0 ) {
        blockIter = blockFirst;
        step = count/2;
        advance(blockIter, step);

        const BtiBlock& block = (*blockIter);
        if ( block.StartPosition <= region.RightPosition ) {
            if ( block.MaxEndPosition > region.LeftPosition ) {
                offset = block.StartOffset;
                break;
            }
            blockFirst = ++blockIter;
            count -= step+1;
        }
        else count = step;
    }

    // if we didn't search "off the end" of the blocks
    if ( blockIter != blockLast ) {

        // "walk back" until we've gone too far
        while ( blockIter != blockFirst ) {
            const BtiBlock& currentBlock = (*blockIter);

            --blockIter;
            const BtiBlock& previousBlock = (*blockIter);
            if ( previousBlock.MaxEndPosition <= region.LeftPosition ) {
                offset = currentBlock.StartOffset;
                found = true;
                break;
            }
        }

        // if we walked all the way to first block, just return that and let the reader's
        // region overlap parsing do the rest
        if ( blockIter == blockFirst ) {
            const BtiBlock& block = (*blockIter);
            offset = block.StartOffset;
            found = true;
        }
    }


    // sets to false if blocks container is empty, or if no matching block could be found
    *hasAlignmentsInRegion = found;
}

// returns whether reference has alignments or no
bool BamToolsIndex::HasAlignments(const int& referenceID) const {
    if ( referenceID < 0 || referenceID >= (int)m_indexFileSummary.size() )
        return false;
    const BtiReferenceSummary& refSummary = m_indexFileSummary.at(referenceID);
    return ( refSummary.NumBlocks > 0 );
}

// pre-allocates space for each reference's summary data
void BamToolsIndex::InitializeFileSummary(const int& numReferences) {
    m_indexFileSummary.clear();
    for ( int i = 0; i < numReferences; ++i )
        m_indexFileSummary.push_back( BtiReferenceSummary() );
}

// returns true if the index stream is open
bool BamToolsIndex::IsDeviceOpen(void) const {
    if ( m_resources.Device == 0 )
        return false;
    return m_resources.Device->IsOpen();
}

// attempts to use index data to jump to @region, returns success/fail
// a "successful" jump indicates no error, but not whether this region has data
//   * thus, the method sets a flag to indicate whether there are alignments
//     available after the jump position
bool BamToolsIndex::Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion) {

    // clear flag
    *hasAlignmentsInRegion = false;

    // skip if invalid reader or not open
    if ( m_reader == 0 || !m_reader->IsOpen() ) {
        SetErrorString("BamToolsIndex::Jump", "could not jump: reader is not open");
        return false;
    }

    // make sure left-bound position is valid
    const RefVector& references = m_reader->GetReferenceData();
    if ( region.LeftPosition > references.at(region.LeftRefID).RefLength ) {
        SetErrorString("BamToolsIndex::Jump", "could not create index: invalid region requested");
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

    // return success/failure of seek
    return m_reader->Seek(offset);
}

// loads existing data from file into memory
bool BamToolsIndex::Load(const std::string& filename) {

    try {

        // attempt to open file (read-only)
        OpenFile(filename, IBamIODevice::ReadOnly);

        // load metadata & generate in-memory summary
        LoadHeader();
        LoadFileSummary();

        // return success
        return true;

    } catch ( BamException& e ) {
        m_errorString = e.what();
        return false;
    }
}

void BamToolsIndex::LoadFileSummary(void) {

    // load number of reference sequences
    int numReferences;
    LoadNumReferences(numReferences);

    // initialize file summary data
    InitializeFileSummary(numReferences);

    // load summary for each reference
    BtiFileSummary::iterator summaryIter = m_indexFileSummary.begin();
    BtiFileSummary::iterator summaryEnd  = m_indexFileSummary.end();
    for ( ; summaryIter != summaryEnd; ++summaryIter )
        LoadReferenceSummary(*summaryIter);
}

void BamToolsIndex::LoadHeader(void) {

    // check BTI file metadata
    CheckMagicNumber();
    CheckVersion();

    // use file's BTI block size to set member variable
    const int64_t numBytesRead = m_resources.Device->Read((char*)&m_blockSize, sizeof(m_blockSize));
    if ( m_isBigEndian ) SwapEndian_32(m_blockSize);
    if ( numBytesRead != sizeof(m_blockSize) )
        throw BamException("BamToolsIndex::LoadHeader", "could not read BTI block size");
}

void BamToolsIndex::LoadNumBlocks(int& numBlocks) {
    const int64_t numBytesRead = m_resources.Device->Read((char*)&numBlocks, sizeof(numBlocks));
    if ( m_isBigEndian ) SwapEndian_32(numBlocks);
    if ( numBytesRead != sizeof(numBlocks) )
        throw BamException("BamToolsIndex::LoadNumBlocks", "could not read number of BTI blocks");
}

void BamToolsIndex::LoadNumReferences(int& numReferences) {
    const int64_t numBytesRead = m_resources.Device->Read((char*)&numReferences, sizeof(numReferences));
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    if ( numBytesRead != sizeof(numReferences) )
        throw BamException("BamToolsIndex::LoadNumReferences", "could not read number of references");
}

void BamToolsIndex::LoadReferenceSummary(BtiReferenceSummary& refSummary) {

    // load number of blocks
    int numBlocks;
    LoadNumBlocks(numBlocks);

    // store block summary data for this reference
    refSummary.NumBlocks = numBlocks;
    refSummary.FirstBlockFilePosition = Tell();

    // skip reference's blocks
    SkipBlocks(numBlocks);
}

void BamToolsIndex::OpenFile(const std::string& filename, IBamIODevice::OpenMode mode) {

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
        throw BamException("BamToolsIndex::OpenFile", message);
    }
}

void BamToolsIndex::ReadBlock(BtiBlock& block) {

    // read in block data members
    int64_t numBytesRead = 0;
    numBytesRead += m_resources.Device->Read((char*)&block.MaxEndPosition, sizeof(block.MaxEndPosition));
    numBytesRead += m_resources.Device->Read((char*)&block.StartOffset,    sizeof(block.StartOffset));
    numBytesRead += m_resources.Device->Read((char*)&block.StartPosition,  sizeof(block.StartPosition));

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_32(block.MaxEndPosition);
        SwapEndian_64(block.StartOffset);
        SwapEndian_32(block.StartPosition);
    }

    // check block read ok
    const int expectedBytes = sizeof(block.MaxEndPosition) +
                              sizeof(block.StartOffset) +
                              sizeof(block.StartPosition);
    if ( numBytesRead != expectedBytes )
        throw BamException("BamToolsIndex::ReadBlock", "could not read block");
}

void BamToolsIndex::ReadBlocks(const BtiReferenceSummary& refSummary, BtiBlockVector& blocks) {

    // prep blocks container
    blocks.clear();
    blocks.reserve(refSummary.NumBlocks);

    // skip to first block entry
    Seek( refSummary.FirstBlockFilePosition, SEEK_SET );

    // read & store block entries
    BtiBlock block;
    for ( int i = 0; i < refSummary.NumBlocks; ++i ) {
        ReadBlock(block);
        blocks.push_back(block);
    }
}

void BamToolsIndex::ReadReferenceEntry(BtiReferenceEntry& refEntry) {

    // return false if refId not valid index in file summary structure
    if ( refEntry.ID < 0 || refEntry.ID >= (int)m_indexFileSummary.size() )
        throw BamException("BamToolsIndex::ReadReferenceEntry", "invalid reference requested");

    // use index summary to assist reading the reference's BTI blocks
    const BtiReferenceSummary& refSummary = m_indexFileSummary.at(refEntry.ID);
    ReadBlocks(refSummary, refEntry.Blocks);
}

void BamToolsIndex::Seek(const int64_t& position, const int origin) {
    if ( !m_resources.Device->Seek(position, origin) )
        throw BamException("BamToolsIndex::Seek", "could not seek in BAI file");
}

void BamToolsIndex::SkipBlocks(const int& numBlocks) {
    Seek( numBlocks*BamToolsIndex::SIZEOF_BLOCK, SEEK_CUR );
}

int64_t BamToolsIndex::Tell(void) const {
    return m_resources.Device->Tell();
}

void BamToolsIndex::WriteBlock(const BtiBlock& block) {

    // copy entry data
    int32_t maxEndPosition = block.MaxEndPosition;
    int64_t startOffset    = block.StartOffset;
    int32_t startPosition  = block.StartPosition;

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_32(maxEndPosition);
        SwapEndian_64(startOffset);
        SwapEndian_32(startPosition);
    }

    // write the reference index entry
    int64_t numBytesWritten = 0;
    numBytesWritten += m_resources.Device->Write((const char*)&maxEndPosition, sizeof(maxEndPosition));
    numBytesWritten += m_resources.Device->Write((const char*)&startOffset,    sizeof(startOffset));
    numBytesWritten += m_resources.Device->Write((const char*)&startPosition,  sizeof(startPosition));

    // check block written ok
    const int expectedBytes = sizeof(maxEndPosition) +
                              sizeof(startOffset) +
                              sizeof(startPosition);
    if ( numBytesWritten != expectedBytes )
        throw BamException("BamToolsIndex::WriteBlock", "could not write BTI block");
}

void BamToolsIndex::WriteBlocks(const BtiBlockVector& blocks) {
    BtiBlockVector::const_iterator blockIter = blocks.begin();
    BtiBlockVector::const_iterator blockEnd  = blocks.end();
    for ( ; blockIter != blockEnd; ++blockIter )
        WriteBlock(*blockIter);
}

void BamToolsIndex::WriteHeader(void) {

    int64_t numBytesWritten = 0 ;

    // write BTI index format 'magic number'
    numBytesWritten += m_resources.Device->Write(BamToolsIndex::BTI_MAGIC, 4);

    // write BTI index format version
    int32_t currentVersion = (int32_t)m_outputVersion;
    if ( m_isBigEndian ) SwapEndian_32(currentVersion);
    numBytesWritten += m_resources.Device->Write((const char*)&currentVersion, sizeof(currentVersion));

    // write block size
    uint32_t blockSize = m_blockSize;
    if ( m_isBigEndian ) SwapEndian_32(blockSize);
    numBytesWritten += m_resources.Device->Write((const char*)&blockSize, sizeof(blockSize));

    // write number of references
    int32_t numReferences = m_indexFileSummary.size();
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    numBytesWritten += m_resources.Device->Write((const char*)&numReferences, sizeof(numReferences));

    // check header written ok
    const int expectedBytes = 4 +
                              sizeof(currentVersion) +
                              sizeof(blockSize) +
                              sizeof(numReferences);
    if ( numBytesWritten != expectedBytes )
        throw BamException("BamToolsIndex::WriteHeader", "could not write BTI header");
}

void BamToolsIndex::WriteReferenceEntry(const BtiReferenceEntry& refEntry) {

    // write number of blocks this reference
    uint32_t numBlocks = refEntry.Blocks.size();
    if ( m_isBigEndian ) SwapEndian_32(numBlocks);
    const int64_t numBytesWritten = m_resources.Device->Write((const char*)&numBlocks, sizeof(numBlocks));
    if ( numBytesWritten != sizeof(numBlocks) )
        throw BamException("BamToolsIndex::WriteReferenceEntry", "could not write number of blocks");

    // write actual block entries
    WriteBlocks(refEntry.Blocks);
}
