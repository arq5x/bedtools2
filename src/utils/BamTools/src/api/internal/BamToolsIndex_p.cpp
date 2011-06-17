// ***************************************************************************
// BamToolsIndex.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 27 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the BamTools index format (".bti")
// ***************************************************************************

#include <api/BamAlignment.h>
#include <api/internal/BamReader_p.h>
#include <api/internal/BamToolsIndex_p.h>
#include <api/internal/BgzfStream_p.h>
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

// static BamToolsIndex constants
const int BamToolsIndex::DEFAULT_BLOCK_LENGTH = 1000;
const string BamToolsIndex::BTI_EXTENSION     = ".bti";
const char* const BamToolsIndex::BTI_MAGIC    = "BTI\1";
const int BamToolsIndex::SIZEOF_BLOCK         = sizeof(int32_t)*2 + sizeof(int64_t);

// ctor
BamToolsIndex::BamToolsIndex(Internal::BamReaderPrivate* reader)
    : BamIndex(reader)
    , m_indexStream(0)
    , m_cacheMode(BamIndex::LimitedIndexCaching)
    , m_blockSize(BamToolsIndex::DEFAULT_BLOCK_LENGTH)
    , m_inputVersion(0)
    , m_outputVersion(BTI_1_2) // latest version - used for writing new index files
{
    m_isBigEndian = BamTools::SystemIsBigEndian();
}

// dtor
BamToolsIndex::~BamToolsIndex(void) {
    CloseFile();
}

bool BamToolsIndex::CheckMagicNumber(void) {

    // check 'magic number' to see if file is BTI index
    char magic[4];
    size_t elementsRead = fread(magic, sizeof(char), 4, m_indexStream);
    if ( elementsRead != 4 ) {
        cerr << "BamToolsIndex ERROR: could not read format 'magic' number" << endl;
        return false;
    }

    if ( strncmp(magic, BamToolsIndex::BTI_MAGIC, 4) != 0 ) {
        cerr << "BamToolsIndex ERROR: invalid format" << endl;
        return false;
    }

    // otherwise ok
    return true;
}

// check index file version, return true if OK
bool BamToolsIndex::CheckVersion(void) {

    // read version from file
    size_t elementsRead = fread(&m_inputVersion, sizeof(m_inputVersion), 1, m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(m_inputVersion);

    // if version is negative, or zero
    if ( m_inputVersion <= 0 ) {
        cerr << "BamToolsIndex ERROR: could not load index file: invalid version."
             << endl;
        return false;
    }

    // if version is newer than can be supported by this version of bamtools
    else if ( m_inputVersion > m_outputVersion ) {
        cerr << "BamToolsIndex ERROR: could not load index file. This version of BamTools does not recognize new index file version"
             << endl
             << "Please update BamTools to a more recent version to support this index file."
             << endl;
        return false;
    }

    // ------------------------------------------------------------------
    // check for deprecated, unsupported versions
    // (typically whose format did not accomodate a particular bug fix)

    else if ( (Version)m_inputVersion == BamToolsIndex::BTI_1_0 ) {
        cerr << "BamToolsIndex ERROR: could not load index file. This version of the index contains a bug related to accessing data near reference ends."
             << endl << endl
             << "Please run 'bamtools index -bti -in yourData.bam' to generate an up-to-date, fixed BTI file."
             << endl << endl;
        return false;
    }

    else if ( (Version)m_inputVersion == BamToolsIndex::BTI_1_1 ) {
        cerr << "BamToolsIndex ERROR: could not load index file. This version of the index contains a bug related to handling empty references."
             << endl << endl
             << "Please run 'bamtools index -bti -in yourData.bam' to generate an up-to-date, fixed BTI file."
             << endl << endl;
        return false;
    }

    // otherwise ok
    else return true;
}

void BamToolsIndex::ClearReferenceEntry(BtiReferenceEntry& refEntry) {
    refEntry.ID = -1;
    refEntry.Blocks.clear();
}

void BamToolsIndex::CloseFile(void) {
    if ( IsFileOpen() )
        fclose(m_indexStream);
    m_indexFileSummary.clear();
}

// builds index from associated BAM file & writes out to index file
bool BamToolsIndex::Create(void) {

    // return false if BamReader is invalid or not open
    if ( m_reader == 0 || !m_reader->IsOpen() ) {
        cerr << "BamToolsIndex ERROR: BamReader is not open"
             << ", aborting index creation" << endl;
        return false;
    }

    // rewind BamReader
    if ( !m_reader->Rewind() ) {
        cerr << "BamToolsIndex ERROR: could not rewind BamReader to create index"
             << ", aborting index creation" << endl;
        return false;
    }

    // open new index file (read & write)
    string indexFilename = m_reader->Filename() + Extension();
    if ( !OpenFile(indexFilename, "w+b") ) {
        cerr << "BamToolsIndex ERROR: could not open ouput index file " << indexFilename
             << ", aborting index creation" << endl;
        return false;
    }

    // initialize BtiFileSummary with number of references
    const int& numReferences = m_reader->GetReferenceCount();
    InitializeFileSummary(numReferences);

    // initialize output file
    bool createdOk = true;
    createdOk &= WriteHeader();

    // index building markers
    int32_t currentBlockCount      = 0;
    int64_t currentAlignmentOffset = m_reader->Tell();
    int32_t blockRefId             = -1;
    int32_t blockMaxEndPosition    = -1;
    int64_t blockStartOffset       = currentAlignmentOffset;
    int32_t blockStartPosition     = -1;

    // plow through alignments, storing index entries
    BamAlignment al;
    BtiReferenceEntry refEntry;
    while ( m_reader->LoadNextAlignment(al) ) {

        // if moved to new reference
        if ( al.RefID != blockRefId ) {

            // if first pass, check:
            if ( currentBlockCount == 0 ) {

                // write any empty references up to (but not including) al.RefID
                for ( int i = 0; i < al.RefID; ++i ) {
                    BtiReferenceEntry emptyEntry(i);
                    createdOk &= WriteReferenceEntry(emptyEntry);
                }
            }

            // not first pass:
            else {

                // store previous BTI block data in reference entry
                BtiBlock block(blockMaxEndPosition, blockStartOffset, blockStartPosition);
                refEntry.Blocks.push_back(block);

                // write reference entry, then clear
                createdOk &= WriteReferenceEntry(refEntry);
                ClearReferenceEntry(refEntry);

                // write any empty references between (but not including) the last blockRefID and current al.RefID
                for ( int i = blockRefId+1; i < al.RefID; ++i ) {
                    BtiReferenceEntry emptyEntry(i);
                    createdOk &= WriteReferenceEntry(emptyEntry);
                }

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
        int32_t alignmentEndPosition = al.GetEndPosition();
        if ( alignmentEndPosition > blockMaxEndPosition )
            blockMaxEndPosition = alignmentEndPosition;

        // if block is full, get offset for next block, reset currentBlockCount
        if ( currentBlockCount == m_blockSize ) {

            // store previous block data in reference entry
            BtiBlock block(blockMaxEndPosition, blockStartOffset, blockStartPosition);
            refEntry.Blocks.push_back(block);

            // update markers
            blockStartOffset  = m_reader->Tell();
            currentBlockCount = 0;
        }

        // not the best name, but for the next iteration, this value will be the offset of the *current* alignment
        // necessary because we won't know if this next alignment is on a new reference until we actually read it
        currentAlignmentOffset = m_reader->Tell();
    }

    // after finishing alignments, if any data was read, check:
    if ( blockRefId >= 0 ) {

        // store last BTI block data in reference entry
        BtiBlock block(blockMaxEndPosition, blockStartOffset, blockStartPosition);
        refEntry.Blocks.push_back(block);

        // write last reference entry, then clear
        createdOk &= WriteReferenceEntry(refEntry);
        ClearReferenceEntry(refEntry);

        // then write any empty references remaining at end of file
        for ( int i = blockRefId+1; i < numReferences; ++i ) {
            BtiReferenceEntry emptyEntry(i);
            createdOk &= WriteReferenceEntry(emptyEntry);
        }
    }

    // rewind reader & return result
    createdOk &= m_reader->Rewind();

    // return result
    return createdOk;
}

// returns format's file extension
const std::string BamToolsIndex::Extension(void) {
    return BamToolsIndex::BTI_EXTENSION;
}

bool BamToolsIndex::GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion) {

    // return false ref ID is not a valid index in file summary data
    if ( region.LeftRefID < 0 || region.LeftRefID >= (int)m_indexFileSummary.size() )
        return false;

    // retrieve reference index data for left bound reference
    BtiReferenceEntry refEntry(region.LeftRefID);
    if ( !ReadReferenceEntry(refEntry) ) {
        cerr << "BamToolsIndex ERROR: could not retrieve index data from BTI file" << endl;
        return false;
    }

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
            if ( block.MaxEndPosition >= region.LeftPosition ) {
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
            if ( previousBlock.MaxEndPosition < region.LeftPosition ) {
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

    // return success
    return true;
}

// returns whether reference has alignments or no
bool BamToolsIndex::HasAlignments(const int& referenceID) const {
    if ( referenceID < 0 || referenceID >= (int)m_indexFileSummary.size() )
        return false;
    const BtiReferenceSummary& refSummary = m_indexFileSummary.at(referenceID);
    return ( refSummary.NumBlocks > 0 );
}

void BamToolsIndex::InitializeFileSummary(const int& numReferences) {
    m_indexFileSummary.clear();
    for ( int i = 0; i < numReferences; ++i )
        m_indexFileSummary.push_back( BtiReferenceSummary() );
}

bool BamToolsIndex::IsFileOpen(void) const {
    return ( m_indexStream != 0 );
}

// attempts to use index data to jump to @region, returns success/fail
// a "successful" jump indicates no error, but not whether this region has data
//   * thus, the method sets a flag to indicate whether there are alignments
//     available after the jump position
bool BamToolsIndex::Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion) {

    // clear flag
    *hasAlignmentsInRegion = false;

    // skip if invalid reader or not open
    if ( m_reader == 0 || !m_reader->IsOpen() )
        return false;

    // make sure left-bound position is valid
    const RefVector& references = m_reader->GetReferenceData();
    if ( region.LeftPosition > references.at(region.LeftRefID).RefLength )
        return false;

    // calculate nearest offset to jump to
    int64_t offset;
    if ( !GetOffset(region, offset, hasAlignmentsInRegion) ) {
        cerr << "BamToolsIndex ERROR: could not jump"
             << ", unable to calculate offset for specified region" << endl;
        return false;
    }

    // return success/failure of seek
    return m_reader->Seek(offset);
}

// loads existing data from file into memory
bool BamToolsIndex::Load(const std::string& filename) {

    // attempt open index file (read-only)
    if ( !OpenFile(filename, "rb") ) {
        cerr << "BamToolsIndex ERROR: could not open input index file " << filename
             << ", aborting index load" << endl;
        return false;
    }

    // attempt to load & validate BTI header data
    if ( !LoadHeader() ) {
        cerr << "BamToolsIndex ERROR: could load header from index file " << filename
             << ", aborting index load" << endl;
        CloseFile();
        return false;
    }

    // attempt to load index file summary
    if ( !LoadFileSummary() ) {
        cerr << "BamToolsIndex ERROR: could not generate a summary of index file " << filename
             << ", aborting index load" << endl;
        CloseFile();
        return false;
    }

    // if we get here, index summary is loaded OK
    return true;
}

bool BamToolsIndex::LoadFileSummary(void) {

    // load number of reference sequences
    int numReferences;
    if ( !LoadNumReferences(numReferences) )
        return false;

    // initialize file summary data
    InitializeFileSummary(numReferences);

    // iterate over reference entries
    bool loadedOk = true;
    BtiFileSummary::iterator summaryIter = m_indexFileSummary.begin();
    BtiFileSummary::iterator summaryEnd  = m_indexFileSummary.end();
    for ( ; summaryIter != summaryEnd; ++summaryIter )
        loadedOk &= LoadReferenceSummary(*summaryIter);

    // return result
    return loadedOk;
}

bool BamToolsIndex::LoadHeader(void) {

    // if invalid format 'magic number'
    if ( !CheckMagicNumber() )
        return false;

    // if invalid BTI version
    if ( !CheckVersion() )
        return false;

    // use file's BTI block size to set member variable
    size_t elementsRead = fread(&m_blockSize, sizeof(m_blockSize), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(m_blockSize);
    return ( elementsRead == 1 );
}

bool BamToolsIndex::LoadNumBlocks(int& numBlocks) {
    size_t elementsRead = 0;
    elementsRead += fread(&numBlocks, sizeof(numBlocks), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numBlocks);
    return ( elementsRead == 1 );
}

bool BamToolsIndex::LoadNumReferences(int& numReferences) {
    size_t elementsRead = 0;
    elementsRead += fread(&numReferences, sizeof(numReferences), 1, m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    return ( elementsRead == 1 );
}

bool BamToolsIndex::LoadReferenceSummary(BtiReferenceSummary& refSummary) {

    // load number of blocks
    int numBlocks;
    if ( !LoadNumBlocks(numBlocks) )
        return false;

    // store block summary data for this reference
    refSummary.NumBlocks = numBlocks;
    refSummary.FirstBlockFilePosition = Tell();

    // skip blocks in index file (and return status)
    return SkipBlocks(numBlocks);
}

bool BamToolsIndex::OpenFile(const std::string& filename, const char* mode) {

    // make sure any previous index file is closed
    CloseFile();

    // attempt to open file
    m_indexStream = fopen(filename.c_str(), mode);
    return IsFileOpen();
}

bool BamToolsIndex::ReadBlock(BtiBlock& block) {

    // read in block data members
    size_t elementsRead = 0;
    elementsRead += fread(&block.MaxEndPosition, sizeof(block.MaxEndPosition), 1, m_indexStream);
    elementsRead += fread(&block.StartOffset,    sizeof(block.StartOffset),    1, m_indexStream);
    elementsRead += fread(&block.StartPosition,  sizeof(block.StartPosition),  1, m_indexStream);

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_32(block.MaxEndPosition);
        SwapEndian_64(block.StartOffset);
        SwapEndian_32(block.StartPosition);
    }

    // return success/failure
    return ( elementsRead == 3 );
}

bool BamToolsIndex::ReadBlocks(const BtiReferenceSummary& refSummary, BtiBlockVector& blocks) {

    // prep blocks container
    blocks.clear();
    blocks.reserve(refSummary.NumBlocks);

    // skip to first block entry
    if ( !Seek( refSummary.FirstBlockFilePosition, SEEK_SET ) ) {
        cerr << "BamToolsIndex ERROR: could not seek to position "
             << refSummary.FirstBlockFilePosition << endl;
        return false;
    }

    // read & store block entries
    bool readOk = true;
    BtiBlock block;
    for ( int i = 0; i < refSummary.NumBlocks; ++i ) {
        readOk &= ReadBlock(block);
        blocks.push_back(block);
    }
    return readOk;
}

bool BamToolsIndex::ReadReferenceEntry(BtiReferenceEntry& refEntry) {

    // return false if refId not valid index in file summary structure
    if ( refEntry.ID < 0 || refEntry.ID >= (int)m_indexFileSummary.size() )
        return false;

    // use index summary to assist reading the reference's BTI blocks
    const BtiReferenceSummary& refSummary = m_indexFileSummary.at(refEntry.ID);
    return ReadBlocks(refSummary, refEntry.Blocks);
}

bool BamToolsIndex::Seek(const int64_t& position, const int& origin) {
    return ( fseek64(m_indexStream, position, origin) == 0 );
}

// change the index caching behavior
void BamToolsIndex::SetCacheMode(const BamIndex::IndexCacheMode& mode) {
    m_cacheMode = mode;
    // do nothing else here ? cache mode will be ignored from now on, most likely
}

bool BamToolsIndex::SkipBlocks(const int& numBlocks) {
    return Seek( numBlocks*BamToolsIndex::SIZEOF_BLOCK, SEEK_CUR );
}

int64_t BamToolsIndex::Tell(void) const {
    return ftell64(m_indexStream);
}

bool BamToolsIndex::WriteBlock(const BtiBlock& block) {

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
    size_t elementsWritten = 0;
    elementsWritten += fwrite(&maxEndPosition, sizeof(maxEndPosition), 1, m_indexStream);
    elementsWritten += fwrite(&startOffset,    sizeof(startOffset),    1, m_indexStream);
    elementsWritten += fwrite(&startPosition,  sizeof(startPosition),  1, m_indexStream);
    return ( elementsWritten == 3 );
}

bool BamToolsIndex::WriteBlocks(const BtiBlockVector& blocks) {
    bool writtenOk = true;
    BtiBlockVector::const_iterator blockIter = blocks.begin();
    BtiBlockVector::const_iterator blockEnd  = blocks.end();
    for ( ; blockIter != blockEnd; ++blockIter )
        writtenOk &= WriteBlock(*blockIter);
    return writtenOk;
}

bool BamToolsIndex::WriteHeader(void) {

    size_t elementsWritten = 0;

    // write BTI index format 'magic number'
    elementsWritten += fwrite(BamToolsIndex::BTI_MAGIC, 1, 4, m_indexStream);

    // write BTI index format version
    int32_t currentVersion = (int32_t)m_outputVersion;
    if ( m_isBigEndian ) SwapEndian_32(currentVersion);
    elementsWritten += fwrite(&currentVersion, sizeof(currentVersion), 1, m_indexStream);

    // write block size
    int32_t blockSize = m_blockSize;
    if ( m_isBigEndian ) SwapEndian_32(blockSize);
    elementsWritten += fwrite(&blockSize, sizeof(blockSize), 1, m_indexStream);

    // write number of references
    int32_t numReferences = m_indexFileSummary.size();
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    elementsWritten += fwrite(&numReferences, sizeof(numReferences), 1, m_indexStream);

    // return success/failure of write
    return ( elementsWritten == 7 );
}

bool BamToolsIndex::WriteReferenceEntry(const BtiReferenceEntry& refEntry) {

    size_t elementsWritten = 0;

    // write number of blocks this reference
    uint32_t numBlocks = refEntry.Blocks.size();
    if ( m_isBigEndian ) SwapEndian_32(numBlocks);
    elementsWritten += fwrite(&numBlocks, sizeof(numBlocks), 1, m_indexStream);

    // write actual block entries
    const bool blocksOk = WriteBlocks(refEntry.Blocks);

    // return success/fail
    return ( elementsWritten == 1) && blocksOk;
}
