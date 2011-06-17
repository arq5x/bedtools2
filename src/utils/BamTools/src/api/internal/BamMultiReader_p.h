// ***************************************************************************
// BamMultiReader_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 13 March 2011 (DB)
// ---------------------------------------------------------------------------
// Functionality for simultaneously reading multiple BAM files
// *************************************************************************

#ifndef BAMMULTIREADER_P_H
#define BAMMULTIREADER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <string>
#include <vector>

namespace BamTools {
namespace Internal {

class IBamMultiMerger;

class BamMultiReaderPrivate {

    // constructor / destructor
    public:
        BamMultiReaderPrivate(void);
        ~BamMultiReaderPrivate(void);

    // public interface
    public:

        // file operations
        void Close(void);
        void CloseFile(const std::string& filename);
        void CloseFiles(const std::vector<std::string>& filenames);
        const std::vector<std::string> Filenames(void) const;
        bool Jump(int refID, int position = 0);
        bool Open(const std::vector<std::string>& filenames);
        bool OpenFile(const std::string& filename);
        void PrintFilenames(void) const;
        bool Rewind(void);
        bool SetRegion(const BamRegion& region);

        // access alignment data
        bool GetNextAlignment(BamAlignment& al);
        bool GetNextAlignmentCore(BamAlignment& al);
        bool HasOpenReaders(void);
        void SetSortOrder(const BamMultiReader::SortOrder& order);

        // access auxiliary data
        SamHeader GetHeader(void) const;
        std::string GetHeaderText(void) const;
        int GetReferenceCount(void) const;
        const BamTools::RefVector GetReferenceData(void) const;
        int GetReferenceID(const std::string& refName) const;

        // BAM index operations
        bool CreateIndexes(const BamIndex::IndexType& type = BamIndex::STANDARD);
        bool HasIndexes(void) const;
        bool LocateIndexes(const BamIndex::IndexType& preferredType = BamIndex::STANDARD);
        bool OpenIndexes(const std::vector<std::string>& indexFilenames);
        void SetIndexCacheMode(const BamIndex::IndexCacheMode mode);

    // 'internal' methods
    public:
        IBamMultiMerger* CreateMergerForCurrentSortOrder(void) const;
        const std::string ExtractReadGroup(const std::string& headerLine) const;
        bool HasAlignmentData(void) const;
        bool LoadNextAlignment(BamAlignment& al);
        BamTools::BamReader* OpenReader(const std::string& filename);
        bool RewindReaders(void);
        void SaveNextAlignment(BamTools::BamReader* reader, BamTools::BamAlignment* alignment);
        const std::vector<std::string> SplitHeaderText(const std::string& headerText) const;
        void UpdateAlignmentCache(void);
        void ValidateReaders(void) const;

    // data members
    public:
        typedef std::pair<BamReader*, BamAlignment*> ReaderAlignment;
        std::vector<ReaderAlignment> m_readers;

        IBamMultiMerger* m_alignments;
        bool m_isCoreMode;
        BamMultiReader::SortOrder m_sortOrder;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMMULTIREADER_P_H
