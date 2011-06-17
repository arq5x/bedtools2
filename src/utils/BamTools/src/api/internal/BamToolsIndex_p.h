// ***************************************************************************
// BamToolsIndex.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 5 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides index operations for the BamTools index format (".bti")
// ***************************************************************************

#ifndef BAMTOOLS_INDEX_FORMAT_H
#define BAMTOOLS_INDEX_FORMAT_H

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
#include <string>
#include <vector>

namespace BamTools {
namespace Internal {

// contains data for each 'block' in a BTI index
struct BtiBlock {

    // data members
    int32_t MaxEndPosition;
    int64_t StartOffset;
    int32_t StartPosition;

    // ctor
    BtiBlock(const int32_t& maxEndPosition = 0,
             const int64_t& startOffset    = 0,
             const int32_t& startPosition  = 0)
        : MaxEndPosition(maxEndPosition)
        , StartOffset(startOffset)
        , StartPosition(startPosition)
    { }
};

// convenience typedef for describing a a list of BTI blocks on a reference
typedef std::vector<BtiBlock> BtiBlockVector;

// contains all fields necessary for building, loading, & writing
// full BTI index data for a single reference
struct BtiReferenceEntry {

    // data members
    int32_t ID;
    BtiBlockVector Blocks;

    // ctor
    BtiReferenceEntry(const int& id = -1)
        : ID(id)
    { }
};

// provides (persistent) summary of BtiReferenceEntry's index data
struct BtiReferenceSummary {

    // data members
    int NumBlocks;
    uint64_t FirstBlockFilePosition;

    // ctor
    BtiReferenceSummary(void)
        : NumBlocks(0)
        , FirstBlockFilePosition(0)
    { }
};

// convenience typedef for describing a full BTI index file summary
typedef std::vector<BtiReferenceSummary> BtiFileSummary;

class BamToolsIndex : public BamIndex {

    // keep a list of any supported versions here
    // (might be useful later to handle any 'legacy' versions if the format changes)
    // listed for example like: BTI_1_0 = 1, BTI_1_1 = 2, BTI_1_2 = 3, BTI_2_0 = 4, and so on
    //
    // so a change introduced in (hypothetical) BTI_1_2 would be handled from then on by:
    //
    // if ( indexVersion >= BTI_1_2 )
    //   do something new
    // else
    //   do the old thing
    enum Version { BTI_1_0 = 1
                 , BTI_1_1
                 , BTI_1_2
                 };

    // ctor & dtor
    public:
        BamToolsIndex(Internal::BamReaderPrivate* reader);
        ~BamToolsIndex(void);

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
        bool CheckVersion(void);
        void CloseFile(void);
        bool IsFileOpen(void) const;
        bool OpenFile(const std::string& filename, const char* mode);
        bool Seek(const int64_t& position, const int& origin);
        int64_t Tell(void) const;

    // internal BTI index building methods
    private:
        void ClearReferenceEntry(BtiReferenceEntry& refEntry);

    // internal random-access methods
    private:
        bool GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion);

    // internal BTI summary data methods
    private:
        void InitializeFileSummary(const int& numReferences);
        bool LoadFileSummary(void);
        bool LoadHeader(void);
        bool LoadNumBlocks(int& numBlocks);
        bool LoadNumReferences(int& numReferences);
        bool LoadReferenceSummary(BtiReferenceSummary& refSummary);
        bool SkipBlocks(const int& numBlocks);

    // internal BTI full index input methods
    private:
        bool ReadBlock(BtiBlock& block);
        bool ReadBlocks(const BtiReferenceSummary& refSummary, BtiBlockVector& blocks);
        bool ReadReferenceEntry(BtiReferenceEntry& refEntry);

    // internal BTI full index output methods
    private:
        bool WriteBlock(const BtiBlock& block);
        bool WriteBlocks(const BtiBlockVector& blocks);
        bool WriteHeader(void);
        bool WriteReferenceEntry(const BtiReferenceEntry& refEntry);

    // data members
    private:
        FILE* m_indexStream;
        bool  m_isBigEndian;
        BamIndex::IndexCacheMode m_cacheMode;
        BtiFileSummary m_indexFileSummary;
        int m_blockSize;
        int32_t m_inputVersion; // Version is serialized as int
        Version m_outputVersion;

    // static constants
    private:
        static const int DEFAULT_BLOCK_LENGTH;
        static const std::string BTI_EXTENSION;
        static const char* const BTI_MAGIC;
        static const int SIZEOF_BLOCK;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMTOOLS_INDEX_FORMAT_H
