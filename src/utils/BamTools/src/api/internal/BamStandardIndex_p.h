// ***************************************************************************
// BamStandardIndex.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the standardized BAM index format (".bai")
// ***************************************************************************

#ifndef BAM_STANDARD_INDEX_FORMAT_H
#define BAM_STANDARD_INDEX_FORMAT_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#include <api/BamAux.h>
#include <api/BamIndex.h>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace BamTools {
namespace Internal {

// -----------------------------------------------------------------------------
// BamStandardIndex data structures

// defines start and end of a contiguous run of alignments
struct BaiAlignmentChunk {

    // data members
    uint64_t Start;
    uint64_t Stop;

    // constructor
    BaiAlignmentChunk(const uint64_t& start = 0,
                      const uint64_t& stop = 0)
        : Start(start)
        , Stop(stop)
    { }
};

// comparison operator (for sorting)
inline
bool operator<(const BaiAlignmentChunk& lhs, const BaiAlignmentChunk& rhs) {
    return lhs.Start < rhs.Start;
}

// convenience typedef for a list of all alignment 'chunks' in a BAI bin
typedef std::vector<BaiAlignmentChunk> BaiAlignmentChunkVector;

// convenience typedef for a map of all BAI bins in a reference (ID => chunks)
typedef std::map<uint32_t, BaiAlignmentChunkVector> BaiBinMap;

// convenience typedef for a list of all 'linear offsets' in a reference
typedef std::vector<uint64_t> BaiLinearOffsetVector;

// contains all fields necessary for building, loading, & writing
// full BAI index data for a single reference
struct BaiReferenceEntry {

    // data members
    int32_t ID;
    BaiBinMap Bins;
    BaiLinearOffsetVector LinearOffsets;

    // ctor
    BaiReferenceEntry(const int32_t& id = -1)
        : ID(id)
    { }
};

// provides (persistent) summary of BaiReferenceEntry's index data
struct BaiReferenceSummary {

    // data members
    int NumBins;
    int NumLinearOffsets;
    uint64_t FirstBinFilePosition;
    uint64_t FirstLinearOffsetFilePosition;

    // ctor
    BaiReferenceSummary(void)
        : NumBins(0)
        , NumLinearOffsets(0)
        , FirstBinFilePosition(0)
        , FirstLinearOffsetFilePosition(0)
    { }
};

// convenience typedef for describing a full BAI index file summary
typedef std::vector<BaiReferenceSummary> BaiFileSummary;

// end BamStandardIndex data structures
// -----------------------------------------------------------------------------

class BamStandardIndex : public BamIndex {

    // ctor & dtor
    public:
        BamStandardIndex(Internal::BamReaderPrivate* reader);
        ~BamStandardIndex(void);

    // BamIndex implementation
    public:
        // builds index from associated BAM file & writes out to index file
        bool Create(void);
        // returns whether reference has alignments or no
        bool HasAlignments(const int& referenceID) const;
        // attempts to use index data to jump to @region, returns success/fail
        // a "successful" jump indicates no error, but not whether this region has data
        //   * thus, the method sets a flag to indicate whether there are alignments
        //     available after the jump position
        bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion);
        // loads existing data from file into memory
        bool Load(const std::string& filename);
        // change the index caching behavior
        void SetCacheMode(const BamIndex::IndexCacheMode& mode);
    public:
        // returns format's file extension
        static const std::string Extension(void);

    // internal file ops
    private:
        bool CheckMagicNumber(void);
        void CloseFile(void);
        bool IsFileOpen(void) const;
        bool OpenFile(const std::string& filename, const char* mode);
        bool Seek(const int64_t& position, const int& origin);
        int64_t Tell(void) const;

    // internal BAI index building methods
    private:
        void ClearReferenceEntry(BaiReferenceEntry& refEntry);
        void SaveAlignmentChunkToBin(BaiBinMap& binMap,
                                     const uint32_t& currentBin,
                                     const uint64_t& currentOffset,
                                     const uint64_t& lastOffset);
        void SaveLinearOffsetEntry(BaiLinearOffsetVector& offsets,
                                   const int& alignmentStartPosition,
                                   const int& alignmentStopPosition,
                                   const uint64_t& lastOffset);

    // internal random-access methods
    private:
        bool AdjustRegion(const BamRegion& region, uint32_t& begin, uint32_t& end);
        void CalculateCandidateBins(const uint32_t& begin,
                                    const uint32_t& end,
                                    std::set<uint16_t>& candidateBins);
        bool CalculateCandidateOffsets(const BaiReferenceSummary& refSummary,
                                       const uint64_t& minOffset,
                                       std::set<uint16_t>& candidateBins,
                                       std::vector<int64_t>& offsets);
        uint64_t CalculateMinOffset(const BaiReferenceSummary& refSummary, const uint32_t& begin);
        bool GetOffsets(const BamRegion& region, std::vector<int64_t>& offsets);
        uint64_t LookupLinearOffset(const BaiReferenceSummary& refSummary, const int& index);

    // internal BAI summary (create/load) methods
    private:
        void ReserveForSummary(const int& numReferences);
        void SaveBinsSummary(const int& refId, const int& numBins);
        void SaveLinearOffsetsSummary(const int& refId, const int& numLinearOffsets);
        bool SkipBins(const int& numBins);
        bool SkipLinearOffsets(const int& numLinearOffsets);
        bool SummarizeBins(BaiReferenceSummary& refSummary);
        bool SummarizeIndexFile(void);
        bool SummarizeLinearOffsets(BaiReferenceSummary& refSummary);
        bool SummarizeReference(BaiReferenceSummary& refSummary);

    // internal BAI full index input methods
    private:
        bool ReadBinID(uint32_t& binId);
        bool ReadBinIntoBuffer(uint32_t& binId, int32_t& numAlignmentChunks);
        bool ReadIntoBuffer(const unsigned int& bytesRequested);
        bool ReadLinearOffset(uint64_t& linearOffset);
        bool ReadNumAlignmentChunks(int& numAlignmentChunks);
        bool ReadNumBins(int& numBins);
        bool ReadNumLinearOffsets(int& numLinearOffsets);
        bool ReadNumReferences(int& numReferences);

    // internal BAI full index output methods
    private:
        void MergeAlignmentChunks(BaiAlignmentChunkVector& chunks);
        void SortLinearOffsets(BaiLinearOffsetVector& linearOffsets);
        bool WriteAlignmentChunk(const BaiAlignmentChunk& chunk);
        bool WriteAlignmentChunks(BaiAlignmentChunkVector& chunks);
        bool WriteBin(const uint32_t& binId, BaiAlignmentChunkVector& chunks);
        bool WriteBins(const int& refId, BaiBinMap& bins);
        bool WriteHeader(void);
        bool WriteLinearOffsets(const int& refId, BaiLinearOffsetVector& linearOffsets);
        bool WriteReferenceEntry(BaiReferenceEntry& refEntry);

    // data members
    private:
        FILE* m_indexStream;
        bool  m_isBigEndian;
        BamIndex::IndexCacheMode m_cacheMode;
        BaiFileSummary m_indexFileSummary;

        // our input buffer
        char* m_buffer;
        unsigned int m_bufferLength;

    // static methods
    private:
        // checks if the buffer is large enough to accomodate the requested size
        static void CheckBufferSize(char*& buffer,
                                    unsigned int& bufferLength,
                                    const unsigned int& requestedBytes);
        // checks if the buffer is large enough to accomodate the requested size
        static void CheckBufferSize(unsigned char*& buffer,
                                    unsigned int& bufferLength,
                                    const unsigned int& requestedBytes);
    // static constants
    private:
        static const int MAX_BIN;
        static const int BAM_LIDX_SHIFT;
        static const std::string BAI_EXTENSION;
        static const char* const BAI_MAGIC;
        static const int SIZEOF_ALIGNMENTCHUNK;
        static const int SIZEOF_BINCORE;
        static const int SIZEOF_LINEAROFFSET;
};

} // namespace Internal
} // namespace BamTools

#endif // BAM_STANDARD_INDEX_FORMAT_H
