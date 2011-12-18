// ***************************************************************************
// BamRandomAccessController_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011(DB)
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

#include "api/BamAux.h"
#include "api/BamIndex.h"

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

    // BamRandomAccessController interface
    public:

        // index methods
        void ClearIndex(void);
        bool CreateIndex(BamReaderPrivate* reader, const BamIndex::IndexType& type);
        bool HasIndex(void) const;
        bool IndexHasAlignmentsForReference(const int& refId);
        bool LocateIndex(BamReaderPrivate* reader, const BamIndex::IndexType& preferredType);
        bool OpenIndex(const std::string& indexFilename, BamReaderPrivate* reader);
        void SetIndex(BamIndex* index);

        // region methods
        void ClearRegion(void);
        bool HasRegion(void) const;
        RegionState AlignmentState(const BamAlignment& alignment) const;
        bool RegionHasAlignments(void) const;
        bool SetRegion(const BamRegion& region, const int& referenceCount);

        // general methods
        void Close(void);
        std::string GetErrorString(void) const;

    // internal methods
    private:
        // adjusts requested region if necessary (depending on where data actually begins)
        void AdjustRegion(const int& referenceCount);
        // error-string handling
        void SetErrorString(const std::string& where, const std::string& what);

    // data members
    private:

        // index data
        BamIndex* m_index;  // owns the index, not a copy - responsible for deleting

        // region data
        BamRegion m_region;
        bool m_hasAlignmentsInRegion;

        // general data
        std::string m_errorString;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMRACONTROLLER_P_H
