// ***************************************************************************
// BamIndex.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 17 August 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index functionality - both for the default (standardized) BAM 
// index format (.bai) as well as a BamTools-specific (nonstandard) index 
// format (.bti).
// ***************************************************************************

#include <cstdio>
#include <cstdlib>
#include <algorithm>
// #include <iostream>
#include <map>
#include "BamIndex.h"
#include "BamReader.h"
#include "BGZF.h"
using namespace std;
using namespace BamTools;

// -------------------------------
// BamIndex implementation

BamIndex::BamIndex(BamTools::BgzfData* bgzf, BamTools::BamReader* reader, bool isBigEndian) 
    : m_BGZF(bgzf)
    , m_reader(reader)
    , m_isBigEndian(isBigEndian)
{ 
    if ( m_reader && m_reader->IsOpen() ) 
        m_references = m_reader->GetReferenceData();
}

bool BamIndex::HasAlignments(const int& referenceID) {
    
    // return false if invalid ID
    if ( (referenceID < 0) || (referenceID >= (int)m_references.size()) ) 
        return false;
    
    // else return status of reference (has alignments?)
    else
        return m_references.at(referenceID).RefHasAlignments;
}

// #########################################################################################
// #########################################################################################

// -------------------------------
// BamDefaultIndex structs & typedefs 
 
namespace BamTools { 

// --------------------------------------------------
// BamDefaultIndex data structures & typedefs
struct Chunk {

    // data members
    uint64_t Start;
    uint64_t Stop;

    // constructor
    Chunk(const uint64_t& start = 0, 
          const uint64_t& stop = 0)
        : Start(start)
        , Stop(stop)
    { }
};

bool ChunkLessThan(const Chunk& lhs, const Chunk& rhs) {
    return lhs.Start < rhs.Start;
}

typedef vector<Chunk> ChunkVector;
typedef map<uint32_t, ChunkVector> BamBinMap;
typedef vector<uint64_t> LinearOffsetVector;

struct ReferenceIndex {
    
    // data members
    BamBinMap Bins;
    LinearOffsetVector Offsets;
    
    // constructor
    ReferenceIndex(const BamBinMap& binMap = BamBinMap(),
                   const LinearOffsetVector& offsets = LinearOffsetVector())
        : Bins(binMap)
        , Offsets(offsets)
    { }
};

typedef vector<ReferenceIndex> BamDefaultIndexData;

} // namespace BamTools
 
// -------------------------------
// BamDefaultIndex implementation
  
struct BamDefaultIndex::BamDefaultIndexPrivate { 
  
    // -------------------------
    // data members
    
    BamDefaultIndexData m_indexData;
    BamDefaultIndex*    m_parent;
    
    // -------------------------
    // ctor & dtor
    
    BamDefaultIndexPrivate(BamDefaultIndex* parent) : m_parent(parent) { }
    ~BamDefaultIndexPrivate(void) { }
    
    // -------------------------
    // internal methods
    
    // calculate bins that overlap region
    int BinsFromRegion(const BamTools::BamRegion& region, const bool isRightBoundSpecified, uint16_t bins[BamTools::MAX_BIN]);
    // saves BAM bin entry for index
    void InsertBinEntry(BamBinMap& binMap, const uint32_t& saveBin, const uint64_t& saveOffset, const uint64_t& lastOffset);
    // saves linear offset entry for index
    void InsertLinearOffset(LinearOffsetVector& offsets, const BamAlignment& bAlignment, const uint64_t& lastOffset);
    // simplifies index by merging 'chunks'
    void MergeChunks(void);
    
};
 
BamDefaultIndex::BamDefaultIndex(BgzfData* bgzf, BamReader* reader, bool isBigEndian)
    : BamIndex(bgzf, reader, isBigEndian)
{
    d = new BamDefaultIndexPrivate(this);
}    

BamDefaultIndex::~BamDefaultIndex(void) {
    d->m_indexData.clear();
    delete d;
    d = 0;
}

// calculate bins that overlap region
int BamDefaultIndex::BamDefaultIndexPrivate::BinsFromRegion(const BamRegion& region, const bool isRightBoundSpecified, uint16_t bins[MAX_BIN]) {
  
    // get region boundaries
    uint32_t begin = (unsigned int)region.LeftPosition;
    uint32_t end;
    
    // if right bound specified AND left&right bounds are on same reference
    // OK to use right bound position
    if ( isRightBoundSpecified && ( region.LeftRefID == region.RightRefID ) )
        end = (unsigned int)region.RightPosition;
    
    // otherwise, use end of left bound reference as cutoff
    else
        end = (unsigned int)m_parent->m_references.at(region.LeftRefID).RefLength - 1;
    
    // initialize list, bin '0' always a valid bin
    int i = 0;
    bins[i++] = 0;

    // get rest of bins that contain this region
    unsigned int k;
    for (k =    1 + (begin>>26); k <=    1 + (end>>26); ++k) { bins[i++] = k; }
    for (k =    9 + (begin>>23); k <=    9 + (end>>23); ++k) { bins[i++] = k; }
    for (k =   73 + (begin>>20); k <=   73 + (end>>20); ++k) { bins[i++] = k; }
    for (k =  585 + (begin>>17); k <=  585 + (end>>17); ++k) { bins[i++] = k; }
    for (k = 4681 + (begin>>14); k <= 4681 + (end>>14); ++k) { bins[i++] = k; }

    // return number of bins stored
    return i;
}

bool BamDefaultIndex::Build(void) { 
  
    // be sure reader & BGZF file are valid & open for reading
    if ( m_reader == 0 || m_BGZF == 0 || !m_BGZF->IsOpen ) 
        return false;

    // move file pointer to beginning of alignments
    m_reader->Rewind();

    // get reference count, reserve index space
    int numReferences = (int)m_references.size();
    for ( int i = 0; i < numReferences; ++i ) {
        d->m_indexData.push_back(ReferenceIndex());
    }
    
    // sets default constant for bin, ID, offset, coordinate variables
    const uint32_t defaultValue = 0xffffffffu;

    // bin data
    uint32_t saveBin(defaultValue);
    uint32_t lastBin(defaultValue);

    // reference ID data
    int32_t saveRefID(defaultValue);
    int32_t lastRefID(defaultValue);

    // offset data
    uint64_t saveOffset = m_BGZF->Tell();
    uint64_t lastOffset = saveOffset;

    // coordinate data
    int32_t lastCoordinate = defaultValue;

    BamAlignment bAlignment;
    while ( m_reader->GetNextAlignmentCore(bAlignment) ) {

        // change of chromosome, save ID, reset bin
        if ( lastRefID != bAlignment.RefID ) {
            lastRefID = bAlignment.RefID;
            lastBin   = defaultValue;
        }

        // if lastCoordinate greater than BAM position - file not sorted properly
        else if ( lastCoordinate > bAlignment.Position ) {
            printf("BAM file not properly sorted:\n");
            printf("Alignment %s : %d > %d on reference (id = %d)", bAlignment.Name.c_str(), lastCoordinate, bAlignment.Position, bAlignment.RefID);
            exit(1);
        }

        // if valid reference && BAM bin spans some minimum cutoff (smaller bin ids span larger regions)
        if ( (bAlignment.RefID >= 0) && (bAlignment.Bin < 4681) ) {

            // save linear offset entry (matched to BAM entry refID)
            ReferenceIndex& refIndex = d->m_indexData.at(bAlignment.RefID);
            LinearOffsetVector& offsets = refIndex.Offsets;
            d->InsertLinearOffset(offsets, bAlignment, lastOffset);
        }

        // if current BamAlignment bin != lastBin, "then possibly write the binning index"
        if ( bAlignment.Bin != lastBin ) {

            // if not first time through
            if ( saveBin != defaultValue ) {

                // save Bam bin entry
                ReferenceIndex& refIndex = d->m_indexData.at(saveRefID);
                BamBinMap& binMap = refIndex.Bins;
                d->InsertBinEntry(binMap, saveBin, saveOffset, lastOffset);
            }

            // update saveOffset
            saveOffset = lastOffset;

            // update bin values
            saveBin = bAlignment.Bin;
            lastBin = bAlignment.Bin;

            // update saveRefID
            saveRefID = bAlignment.RefID;

            // if invalid RefID, break out (why?)
            if ( saveRefID < 0 ) { break; }
        }

        // make sure that current file pointer is beyond lastOffset
        if ( m_BGZF->Tell() <= (int64_t)lastOffset ) {
            printf("Error in BGZF offsets.\n");
            exit(1);
        }

        // update lastOffset
        lastOffset = m_BGZF->Tell();

        // update lastCoordinate
        lastCoordinate = bAlignment.Position;
    }

    // save any leftover BAM data (as long as refID is valid)
    if ( saveRefID >= 0 ) {
        // save Bam bin entry
        ReferenceIndex& refIndex = d->m_indexData.at(saveRefID);
        BamBinMap& binMap = refIndex.Bins;
        d->InsertBinEntry(binMap, saveBin, saveOffset, lastOffset);
    }

    // simplify index by merging chunks
    d->MergeChunks();
    
    // iterate through references in index
    // store whether reference has data &
    // sort offsets in linear offset vector
    BamDefaultIndexData::iterator indexIter = d->m_indexData.begin();
    BamDefaultIndexData::iterator indexEnd  = d->m_indexData.end();
    for ( int i = 0; indexIter != indexEnd; ++indexIter, ++i ) {

        // get reference index data
        ReferenceIndex& refIndex = (*indexIter);
        BamBinMap& binMap = refIndex.Bins;
        LinearOffsetVector& offsets = refIndex.Offsets;

        // store whether reference has alignments or no
        m_references[i].RefHasAlignments = ( binMap.size() > 0 );

        // sort linear offsets
        sort(offsets.begin(), offsets.end());
    }

    // rewind file pointer to beginning of alignments, return success/fail
    return m_reader->Rewind();
}

bool BamDefaultIndex::GetOffsets(const BamRegion& region, const bool isRightBoundSpecified, vector<int64_t>& offsets) { 
  
    // calculate which bins overlap this region
    uint16_t* bins = (uint16_t*)calloc(MAX_BIN, 2);
    int numBins = d->BinsFromRegion(region, isRightBoundSpecified, bins);

    // get bins for this reference
    const ReferenceIndex& refIndex = d->m_indexData.at(region.LeftRefID);
    const BamBinMap& binMap        = refIndex.Bins;

    // get minimum offset to consider
    const LinearOffsetVector& linearOffsets = refIndex.Offsets;
    uint64_t minOffset = ( (unsigned int)(region.LeftPosition>>BAM_LIDX_SHIFT) >= linearOffsets.size() ) ? 0 : linearOffsets.at(region.LeftPosition>>BAM_LIDX_SHIFT);

    // store all alignment 'chunk' starts (file offsets) for bins in this region
    for ( int i = 0; i < numBins; ++i ) {
      
        const uint16_t binKey = bins[i];
        map<uint32_t, ChunkVector>::const_iterator binIter = binMap.find(binKey);
        if ( (binIter != binMap.end()) && ((*binIter).first == binKey) ) {

            const ChunkVector& chunks = (*binIter).second;
            std::vector<Chunk>::const_iterator chunksIter = chunks.begin();
            std::vector<Chunk>::const_iterator chunksEnd  = chunks.end();
            for ( ; chunksIter != chunksEnd; ++chunksIter) {
              
                // if valid chunk found, store its file offset
                const Chunk& chunk = (*chunksIter);
                if ( chunk.Stop > minOffset )
                    offsets.push_back( chunk.Start );
            }
        }
    }

    // clean up memory
    free(bins);

    // sort the offsets before returning
    sort(offsets.begin(), offsets.end());
    
    // return whether any offsets were found
    return ( offsets.size() != 0 );
}

// saves BAM bin entry for index
void BamDefaultIndex::BamDefaultIndexPrivate::InsertBinEntry(BamBinMap&      binMap,
                                     const uint32_t& saveBin,
                                     const uint64_t& saveOffset,
                                     const uint64_t& lastOffset)
{
    // look up saveBin
    BamBinMap::iterator binIter = binMap.find(saveBin);

    // create new chunk
    Chunk newChunk(saveOffset, lastOffset);

    // if entry doesn't exist
    if ( binIter == binMap.end() ) {
        ChunkVector newChunks;
        newChunks.push_back(newChunk);
        binMap.insert( pair<uint32_t, ChunkVector>(saveBin, newChunks));
    }

    // otherwise
    else {
        ChunkVector& binChunks = (*binIter).second;
        binChunks.push_back( newChunk );
    }
}

// saves linear offset entry for index
void BamDefaultIndex::BamDefaultIndexPrivate::InsertLinearOffset(LinearOffsetVector& offsets,
                                        const BamAlignment& bAlignment,
                                        const uint64_t&     lastOffset)
{
    // get converted offsets
    int beginOffset = bAlignment.Position >> BAM_LIDX_SHIFT;
    int endOffset   = (bAlignment.GetEndPosition() - 1) >> BAM_LIDX_SHIFT;

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

bool BamDefaultIndex::Load(const string& filename)  { 
    
    // open index file, abort on error
    FILE* indexStream = fopen(filename.c_str(), "rb");
    if( !indexStream ) {
        printf("ERROR: Unable to open the BAM index file %s for reading.\n", filename.c_str());
        return false;
    }

    // set placeholder to receive input byte count (suppresses compiler warnings)
    size_t elementsRead = 0;
        
    // see if index is valid BAM index
    char magic[4];
    elementsRead = fread(magic, 1, 4, indexStream);
    if ( strncmp(magic, "BAI\1", 4) ) {
        printf("Problem with index file - invalid format.\n");
        fclose(indexStream);
        return false;
    }

    // get number of reference sequences
    uint32_t numRefSeqs;
    elementsRead = fread(&numRefSeqs, 4, 1, indexStream);
    if ( m_isBigEndian ) { SwapEndian_32(numRefSeqs); }
    
    // intialize space for BamDefaultIndexData data structure
    d->m_indexData.reserve(numRefSeqs);

    // iterate over reference sequences
    for ( unsigned int i = 0; i < numRefSeqs; ++i ) {

        // get number of bins for this reference sequence
        int32_t numBins;
        elementsRead = fread(&numBins, 4, 1, indexStream);
        if ( m_isBigEndian ) { SwapEndian_32(numBins); }

        if ( numBins > 0 ) {
            RefData& refEntry = m_references[i];
            refEntry.RefHasAlignments = true;
        }

        // intialize BinVector
        BamBinMap binMap;

        // iterate over bins for that reference sequence
        for ( int j = 0; j < numBins; ++j ) {

            // get binID
            uint32_t binID;
            elementsRead = fread(&binID, 4, 1, indexStream);

            // get number of regionChunks in this bin
            uint32_t numChunks;
            elementsRead = fread(&numChunks, 4, 1, indexStream);

            if ( m_isBigEndian ) { 
              SwapEndian_32(binID);
              SwapEndian_32(numChunks);
            }
            
            // intialize ChunkVector
            ChunkVector regionChunks;
            regionChunks.reserve(numChunks);

            // iterate over regionChunks in this bin
            for ( unsigned int k = 0; k < numChunks; ++k ) {

                // get chunk boundaries (left, right)
                uint64_t left;
                uint64_t right;
                elementsRead = fread(&left, 8, 1, indexStream);
                elementsRead = fread(&right, 8, 1, indexStream);

                if ( m_isBigEndian ) {
                    SwapEndian_64(left);
                    SwapEndian_64(right);
                }
                
                // save ChunkPair
                regionChunks.push_back( Chunk(left, right) );
            }

            // sort chunks for this bin
            sort( regionChunks.begin(), regionChunks.end(), ChunkLessThan );

            // save binID, chunkVector for this bin
            binMap.insert( pair<uint32_t, ChunkVector>(binID, regionChunks) );
        }

        // load linear index for this reference sequence

        // get number of linear offsets
        int32_t numLinearOffsets;
        elementsRead = fread(&numLinearOffsets, 4, 1, indexStream);
        if ( m_isBigEndian ) { SwapEndian_32(numLinearOffsets); }

        // intialize LinearOffsetVector
        LinearOffsetVector offsets;
        offsets.reserve(numLinearOffsets);

        // iterate over linear offsets for this reference sequeence
        uint64_t linearOffset;
        for ( int j = 0; j < numLinearOffsets; ++j ) {
            // read a linear offset & store
            elementsRead = fread(&linearOffset, 8, 1, indexStream);
            if ( m_isBigEndian ) { SwapEndian_64(linearOffset); }
            offsets.push_back(linearOffset);
        }

        // sort linear offsets
        sort( offsets.begin(), offsets.end() );

        // store index data for that reference sequence
        d->m_indexData.push_back( ReferenceIndex(binMap, offsets) );
    }

    // close index file (.bai) and return
    fclose(indexStream);
    return true;
}

// merges 'alignment chunks' in BAM bin (used for index building)
void BamDefaultIndex::BamDefaultIndexPrivate::MergeChunks(void) {

    // iterate over reference enties
    BamDefaultIndexData::iterator indexIter = m_indexData.begin();
    BamDefaultIndexData::iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {

        // get BAM bin map for this reference
        ReferenceIndex& refIndex = (*indexIter);
        BamBinMap& bamBinMap = refIndex.Bins;

        // iterate over BAM bins
        BamBinMap::iterator binIter = bamBinMap.begin();
        BamBinMap::iterator binEnd  = bamBinMap.end();
        for ( ; binIter != binEnd; ++binIter ) {

            // get chunk vector for this bin
            ChunkVector& binChunks = (*binIter).second;
            if ( binChunks.size() == 0 ) { continue; }

            ChunkVector mergedChunks;
            mergedChunks.push_back( binChunks[0] );

            // iterate over chunks
            int i = 0;
            ChunkVector::iterator chunkIter = binChunks.begin();
            ChunkVector::iterator chunkEnd  = binChunks.end();
            for ( ++chunkIter; chunkIter != chunkEnd; ++chunkIter) {

                // get 'currentChunk' based on numeric index
                Chunk& currentChunk = mergedChunks[i];

                // get iteratorChunk based on vector iterator
                Chunk& iteratorChunk = (*chunkIter);

                // if currentChunk.Stop(shifted) == iterator Chunk.Start(shifted)
                if ( currentChunk.Stop>>16 == iteratorChunk.Start>>16 ) {

                    // set currentChunk.Stop to iteratorChunk.Stop
                    currentChunk.Stop = iteratorChunk.Stop;
                }

                // otherwise
                else {
                    // set currentChunk + 1 to iteratorChunk
                    mergedChunks.push_back(iteratorChunk);
                    ++i;
                }
            }

            // saved merged chunk vector
            (*binIter).second = mergedChunks;
        }
    }
}

// writes in-memory index data out to file 
// N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
bool BamDefaultIndex::Write(const std::string& bamFilename) { 

    string indexFilename = bamFilename + ".bai";
    FILE* indexStream = fopen(indexFilename.c_str(), "wb");
    if ( indexStream == 0 ) {
        printf("ERROR: Could not open file to save index.\n");
        return false;
    }

    // write BAM index header
    fwrite("BAI\1", 1, 4, indexStream);

    // write number of reference sequences
    int32_t numReferenceSeqs = d->m_indexData.size();
    if ( m_isBigEndian ) { SwapEndian_32(numReferenceSeqs); }
    fwrite(&numReferenceSeqs, 4, 1, indexStream);

    // iterate over reference sequences
    BamDefaultIndexData::const_iterator indexIter = d->m_indexData.begin();
    BamDefaultIndexData::const_iterator indexEnd  = d->m_indexData.end();
    for ( ; indexIter != indexEnd; ++ indexIter ) {

        // get reference index data
        const ReferenceIndex& refIndex = (*indexIter);
        const BamBinMap& binMap = refIndex.Bins;
        const LinearOffsetVector& offsets = refIndex.Offsets;

        // write number of bins
        int32_t binCount = binMap.size();
        if ( m_isBigEndian ) { SwapEndian_32(binCount); }
        fwrite(&binCount, 4, 1, indexStream);

        // iterate over bins
        BamBinMap::const_iterator binIter = binMap.begin();
        BamBinMap::const_iterator binEnd  = binMap.end();
        for ( ; binIter != binEnd; ++binIter ) {

            // get bin data (key and chunk vector)
            uint32_t binKey = (*binIter).first;
            const ChunkVector& binChunks = (*binIter).second;

            // save BAM bin key
            if ( m_isBigEndian ) { SwapEndian_32(binKey); }
            fwrite(&binKey, 4, 1, indexStream);

            // save chunk count
            int32_t chunkCount = binChunks.size();
            if ( m_isBigEndian ) { SwapEndian_32(chunkCount); }
            fwrite(&chunkCount, 4, 1, indexStream);

            // iterate over chunks
            ChunkVector::const_iterator chunkIter = binChunks.begin();
            ChunkVector::const_iterator chunkEnd  = binChunks.end();
            for ( ; chunkIter != chunkEnd; ++chunkIter ) {

                // get current chunk data
                const Chunk& chunk    = (*chunkIter);
                uint64_t start = chunk.Start;
                uint64_t stop  = chunk.Stop;

                if ( m_isBigEndian ) {
                    SwapEndian_64(start);
                    SwapEndian_64(stop);
                }
                
                // save chunk offsets
                fwrite(&start, 8, 1, indexStream);
                fwrite(&stop,  8, 1, indexStream);
            }
        }

        // write linear offsets size
        int32_t offsetSize = offsets.size();
        if ( m_isBigEndian ) { SwapEndian_32(offsetSize); }
        fwrite(&offsetSize, 4, 1, indexStream);

        // iterate over linear offsets
        LinearOffsetVector::const_iterator offsetIter = offsets.begin();
        LinearOffsetVector::const_iterator offsetEnd  = offsets.end();
        for ( ; offsetIter != offsetEnd; ++offsetIter ) {

            // write linear offset value
            uint64_t linearOffset = (*offsetIter);
            if ( m_isBigEndian ) { SwapEndian_64(linearOffset); }
            fwrite(&linearOffset, 8, 1, indexStream);
        }
    }

    // flush buffer, close file, and return success
    fflush(indexStream);
    fclose(indexStream);
    return true;
}

// #########################################################################################
// #########################################################################################

// -------------------------------------
// BamToolsIndex implementation

namespace BamTools {
  
struct BamToolsIndexEntry {
    
    // data members
    int64_t Offset;
    int RefID;
    int Position;
    
    // ctor
    BamToolsIndexEntry(const uint64_t& offset = 0,
                       const int& id = -1,
                       const int& position = -1)
        : Offset(offset)
        , RefID(id)
        , Position(position)
    { }
};

typedef vector<BamToolsIndexEntry> BamToolsIndexData;
  
} // namespace BamTools

struct BamToolsIndex::BamToolsIndexPrivate {
  
    // -------------------------
    // data members
    BamToolsIndexData m_indexData;
    BamToolsIndex*    m_parent;
    int32_t           m_blockSize;
    
    // -------------------------
    // ctor & dtor
    
    BamToolsIndexPrivate(BamToolsIndex* parent) 
        : m_parent(parent)
        , m_blockSize(1000)
    { }
    
    ~BamToolsIndexPrivate(void) { }
    
    // -------------------------
    // internal methods
};

BamToolsIndex::BamToolsIndex(BgzfData* bgzf, BamReader* reader, bool isBigEndian)
    : BamIndex(bgzf, reader, isBigEndian)
{ 
    d = new BamToolsIndexPrivate(this);
}    

BamToolsIndex::~BamToolsIndex(void) { 
    delete d;
    d = 0;
}

bool BamToolsIndex::Build(void) { 
  
    // be sure reader & BGZF file are valid & open for reading
    if ( m_reader == 0 || m_BGZF == 0 || !m_BGZF->IsOpen ) 
        return false;

    // move file pointer to beginning of alignments
    m_reader->Rewind();
    
    // plow through alignments, store block offsets
    int32_t currentBlockCount  = 0;
    int64_t blockStartOffset   = m_BGZF->Tell();
    int     blockStartId       = -1;
    int     blockStartPosition = -1;
    BamAlignment al;
    while ( m_reader->GetNextAlignmentCore(al) ) {
        
        // set reference flag
        m_references[al.RefID].RefHasAlignments = true;
      
        // if beginning of block, save first alignment's refID & position
        if ( currentBlockCount == 0 ) {
            blockStartId = al.RefID;
            blockStartPosition = al.Position;
        }
      
        // increment block counter
        ++currentBlockCount;
        
        // if block is full, get offset for next block, reset currentBlockCount
        if ( currentBlockCount == d->m_blockSize ) {
          
            d->m_indexData.push_back( BamToolsIndexEntry(blockStartOffset, blockStartId, blockStartPosition) );
            blockStartOffset = m_BGZF->Tell();
            currentBlockCount = 0;
        }
    }
    
    return m_reader->Rewind();
}

// N.B. - ignores isRightBoundSpecified
bool BamToolsIndex::GetOffsets(const BamRegion& region, const bool isRightBoundSpecified, vector<int64_t>& offsets) { 
  
    // return false if no index data present 
    if ( d->m_indexData.empty() ) return false;
  
    // clear any prior data
    offsets.clear();
  
    // calculate nearest index to jump to
    int64_t previousOffset = -1;
    BamToolsIndexData::const_iterator indexIter = d->m_indexData.begin();
    BamToolsIndexData::const_iterator indexEnd  = d->m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {
     
        const BamToolsIndexEntry& entry = (*indexIter);
        
        // check if we are 'past' beginning of desired region
        // if so, we will break out & use previously stored offset
        if ( entry.RefID > region.LeftRefID ) break;
        if ( (entry.RefID == region.LeftRefID) && (entry.Position > region.LeftPosition) ) break;
        
        // not past desired region, so store current entry offset in previousOffset
        previousOffset = entry.Offset;
    }
  
    // no index was found
    if ( previousOffset == -1 ) 
        return false;
    
    // store offset & return success
    offsets.push_back(previousOffset);
    return true; 
}

bool BamToolsIndex::Load(const string& filename) { 
  
    // open index file, abort on error
    FILE* indexStream = fopen(filename.c_str(), "rb");
    if( !indexStream ) {
        printf("ERROR: Unable to open the BAM index file %s for reading.\n", filename.c_str());
        return false;
    }

    // set placeholder to receive input byte count (suppresses compiler warnings)
    size_t elementsRead = 0;
        
    // see if index is valid BAM index
    char magic[4];
    elementsRead = fread(magic, 1, 4, indexStream);
    if ( strncmp(magic, "BTI\1", 4) ) {
        printf("Problem with index file - invalid format.\n");
        fclose(indexStream);
        return false;
    }

    // read in block size
    elementsRead = fread(&d->m_blockSize, sizeof(d->m_blockSize), 1, indexStream);
    if ( m_isBigEndian ) { SwapEndian_32(d->m_blockSize); }
    
    // read in number of offsets
    uint32_t numOffsets;
    elementsRead = fread(&numOffsets, sizeof(numOffsets), 1, indexStream);
    if ( m_isBigEndian ) { SwapEndian_32(numOffsets); }
    
    // reserve space for index data
    d->m_indexData.reserve(numOffsets);

    // iterate over index entries
    for ( unsigned int i = 0; i < numOffsets; ++i ) {
      
        uint64_t offset;
        int id;
        int position;
        
        // read in data
        elementsRead = fread(&offset, sizeof(offset), 1, indexStream);
        elementsRead = fread(&id, sizeof(id), 1, indexStream);
        elementsRead = fread(&position, sizeof(position), 1, indexStream);
        
        // swap endian-ness if necessary
        if ( m_isBigEndian ) {
            SwapEndian_64(offset);
            SwapEndian_32(id);
            SwapEndian_32(position);
        }
        
        // save reference index entry
        d->m_indexData.push_back( BamToolsIndexEntry(offset, id, position) );
        
        // set reference flag
        m_references[id].RefHasAlignments = true;       // what about sparse references? wont be able to set flag?
    }

    // close index file and return
    fclose(indexStream);
    return true;
}

// writes in-memory index data out to file 
// N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
bool BamToolsIndex::Write(const std::string& bamFilename) { 
    
    string indexFilename = bamFilename + ".bti";
    FILE* indexStream = fopen(indexFilename.c_str(), "wb");
    if ( indexStream == 0 ) {
        printf("ERROR: Could not open file to save index.\n");
        return false;
    }

    // write BAM index header
    fwrite("BTI\1", 1, 4, indexStream);

    // write block size
    int32_t blockSize = d->m_blockSize;
    if ( m_isBigEndian ) { SwapEndian_32(blockSize); }
    fwrite(&blockSize, sizeof(blockSize), 1, indexStream);
    
    // write number of offset entries
    uint32_t numOffsets = d->m_indexData.size();
    if ( m_isBigEndian ) { SwapEndian_32(numOffsets); }
    fwrite(&numOffsets, sizeof(numOffsets), 1, indexStream);
    
    // iterate over offset entries
    BamToolsIndexData::const_iterator indexIter = d->m_indexData.begin();
    BamToolsIndexData::const_iterator indexEnd  = d->m_indexData.end();
    for ( ; indexIter != indexEnd; ++ indexIter ) {

        // get reference index data
        const BamToolsIndexEntry& entry = (*indexIter);
        
        // copy entry data
        uint64_t offset = entry.Offset;
        int id = entry.RefID;
        int position = entry.Position;
        
        // swap endian-ness if necessary
        if ( m_isBigEndian ) {
            SwapEndian_64(offset);
            SwapEndian_32(id);
            SwapEndian_32(position);
        }
        
        // write the reference index entry
        fwrite(&offset,   sizeof(offset), 1, indexStream);
        fwrite(&id,       sizeof(id), 1, indexStream);
        fwrite(&position, sizeof(position), 1, indexStream);
    }

    // flush file buffer, close file, and return success
    fflush(indexStream);
    fclose(indexStream);
    return true;
}
