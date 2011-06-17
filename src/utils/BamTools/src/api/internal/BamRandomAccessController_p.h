// ***************************************************************************
// BamRandomAccessController_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 24 February 2011(DB)
// ---------------------------------------------------------------------------
// Manages random access operations in a BAM file
// ***************************************************************************

#ifndef BAMRACONTROLLER_P_H
#define BAMRACONTROLLER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include <api/BamAux.h>
#include <api/BamIndex.h>

namespace BamTools {

class BamAlignment;

namespace Internal {

class BamReaderPrivate;

class BamRandomAccessController {

    // enums
    public: enum RegionState { BeforeRegion = 0
                             , OverlapsRegion
                             , AfterRegion
                             };

    // ctor & dtor
    public:
        BamRandomAccessController(void);
        ~BamRandomAccessController(void);

    // general interface
    public:
        void Close(void);

    // index operations
    public:
        //
        void ClearIndex(void);
        bool CreateIndex(BamReaderPrivate* reader, const BamIndex::IndexType& type);
        bool HasIndex(void) const;
        bool IndexHasAlignmentsForReference(const int& refId);
        bool LocateIndex(BamReaderPrivate* reader, const BamIndex::IndexType& preferredType);
        bool OpenIndex(const std::string& indexFilename, BamReaderPrivate* reader);
        void SetIndex(BamIndex* index);
        void SetIndexCacheMode(const BamIndex::IndexCacheMode& mode);

    // region operations
    public:
        void ClearRegion(void);
        bool HasRegion(void) const;
        RegionState AlignmentState(const BamAlignment& alignment) const;
        bool RegionHasAlignments(void) const;
        bool SetRegion(BamReaderPrivate* reader,
                       const BamRegion& region,
                       const int& referenceCount);

    // 'internal' methods
    public:
        // adjusts requested region if necessary (depending on where data actually begins)
        void AdjustRegion(const int& referenceCount);

    // data members
    private:

        // index data
        BamIndex* m_index;  // owns index, not a copy - responsible for deleting
        BamIndex::IndexCacheMode m_indexCacheMode;

        // region data
        BamRegion m_region;
        bool m_hasAlignmentsInRegion;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMRACONTROLLER_P_H
