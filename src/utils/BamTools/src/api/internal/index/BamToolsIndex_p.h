// ***************************************************************************
// BamToolsIndex.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 November 2011 (DB)
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

#include "api/BamAux.h"
#include "api/BamIndex.h"
#include "api/IBamIODevice.h"
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
    // so a change introduced in BTI_1_2 may be handled from then on by:
    //
    // if ( indexVersion >= BTI_1_2 )
    //   do something new
    // else
    //   do the old thing
    enum Version { BTI_1_0 = 1
                 , BTI_1_1
                 , BTI_1_2
                 , BTI_2_0
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
        BamIndex::IndexType Type(void) const { return BamIndex::BAMTOOLS; }
    public:
        // returns format's file extension
        static const std::string Extension(void);

    // internal methods
    private:

        // index file ops
        void CheckMagicNumber(void);
        void CheckVersion(void);
        void CloseFile(void);
        bool IsDeviceOpen(void) const;
        void OpenFile(const std::string& filename, IBamIODevice::OpenMode mode);
        void Seek(const int64_t& position, const int origin);
        int64_t Tell(void) const;

        // index-creation methods
        void ClearReferenceEntry(BtiReferenceEntry& refEntry);
        void WriteBlock(const BtiBlock& block);
        void WriteBlocks(const BtiBlockVector& blocks);
        void WriteHeader(void);
        void WriteReferenceEntry(const BtiReferenceEntry& refEntry);

        // random-access methods
        void GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion);
        void ReadBlock(BtiBlock& block);
        void ReadBlocks(const BtiReferenceSummary& refSummary, BtiBlockVector& blocks);
        void ReadReferenceEntry(BtiReferenceEntry& refEntry);

        // BTI summary data methods
        void InitializeFileSummary(const int& numReferences);
        void LoadFileSummary(void);
        void LoadHeader(void);
        void LoadNumBlocks(int& numBlocks);
        void LoadNumReferences(int& numReferences);
        void LoadReferenceSummary(BtiReferenceSummary& refSummary);
        void SkipBlocks(const int& numBlocks);

    // data members
    private:
        bool  m_isBigEndian;
        BtiFileSummary m_indexFileSummary;
        uint32_t m_blockSize;
        int32_t m_inputVersion; // Version is serialized as int
        Version m_outputVersion;

        struct RaiiWrapper {
            IBamIODevice* Device;
            RaiiWrapper(void);
            ~RaiiWrapper(void);
        };
        RaiiWrapper m_resources;

    // static constants
    private:
        static const uint32_t DEFAULT_BLOCK_LENGTH;
        static const std::string BTI_EXTENSION;
        static const char* const BTI_MAGIC;
        static const int SIZEOF_BLOCK;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMTOOLS_INDEX_FORMAT_H
