// ***************************************************************************
// BamRandomAccessController_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011(DB)
// ---------------------------------------------------------------------------
// Manages random access operations in a BAM file
// **************************************************************************

#include "api/BamIndex.h"
#include "api/internal/bam/BamRandomAccessController_p.h"
#include "api/internal/bam/BamReader_p.h"
#include "api/internal/index/BamIndexFactory_p.h"
#include "api/internal/utils/BamException_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cassert>
#include <sstream>
using namespace std;

BamRandomAccessController::BamRandomAccessController(void)
    : m_index(0)
    , m_hasAlignmentsInRegion(true)
{ }

BamRandomAccessController::~BamRandomAccessController(void) {
    Close();
}

void BamRandomAccessController::AdjustRegion(const int& referenceCount) {

    // skip if no index available
    if ( m_index == 0 )
        return;

    // see if any references in region have alignments
    m_hasAlignmentsInRegion = false;
    int currentId = m_region.LeftRefID;
    const int rightBoundRefId = ( m_region.isRightBoundSpecified() ? m_region.RightRefID : referenceCount - 1 );
    while ( currentId <= rightBoundRefId ) {
        m_hasAlignmentsInRegion = m_index->HasAlignments(currentId);
        if ( m_hasAlignmentsInRegion ) break;
        ++currentId;
    }

    // if no data found on any reference in region
    if ( !m_hasAlignmentsInRegion )
        return;

    // if left bound of desired region had no data, use first reference that had data
    // otherwise, leave requested region as-is
    if ( currentId != m_region.LeftRefID ) {
        m_region.LeftRefID = currentId;
        m_region.LeftPosition = 0;
    }
}

// returns alignments' "RegionState": { Before|Overlaps|After } current region
BamRandomAccessController::RegionState
BamRandomAccessController::AlignmentState(const BamAlignment& alignment) const {

    // if region has no left bound at all
    if ( !m_region.isLeftBoundSpecified() )
        return OverlapsRegion;

    // handle unmapped reads - return AFTER region to halt processing
    if ( alignment.RefID == -1 )
        return AfterRegion;

    // if alignment is on any reference before left bound reference
    if ( alignment.RefID < m_region.LeftRefID )
        return BeforeRegion;

    // if alignment is on left bound reference
    else if ( alignment.RefID == m_region.LeftRefID ) {

        // if alignment starts at or after left bound position
        if ( alignment.Position >= m_region.LeftPosition) {

            if ( m_region.isRightBoundSpecified() &&             // right bound is specified AND
                 m_region.LeftRefID == m_region.RightRefID &&    // left & right bounds on same reference AND
                 alignment.Position >= m_region.RightPosition )  // alignment starts on or after right bound position
                return AfterRegion;

            // otherwise, alignment overlaps region
            else return OverlapsRegion;
        }

        // alignment starts before left bound position
        else {

            // if alignment overlaps left bound position
            if ( alignment.GetEndPosition() > m_region.LeftPosition )
                return OverlapsRegion;
            else
                return BeforeRegion;
        }
    }

    // otherwise alignment is on a reference after left bound reference
    else {

        // if region has a right bound
        if ( m_region.isRightBoundSpecified() ) {

            // alignment is on any reference between boundaries
            if ( alignment.RefID < m_region.RightRefID )
                return OverlapsRegion;

            // alignment is on any reference after right boundary
            else if ( alignment.RefID > m_region.RightRefID )
                return AfterRegion;

            // alignment is on right bound reference
            else {

                // if alignment starts before right bound position
                if ( alignment.Position < m_region.RightPosition )
                    return OverlapsRegion;
                else
                    return AfterRegion;
            }
        }

        // otherwise, alignment starts after left bound and there is no right bound given
        else return OverlapsRegion;
    }
}

void BamRandomAccessController::Close(void) {
    ClearIndex();
    ClearRegion();
}

void BamRandomAccessController::ClearIndex(void) {
    if ( m_index ) {
        delete m_index;
        m_index = 0;
    }
}

void BamRandomAccessController::ClearRegion(void) {
    m_region.clear();
    m_hasAlignmentsInRegion = true;
}

bool BamRandomAccessController::CreateIndex(BamReaderPrivate* reader,
                                            const BamIndex::IndexType& type)
{
    // skip if reader is invalid
    assert(reader);
    if ( !reader->IsOpen() ) {
        SetErrorString("BamRandomAccessController::CreateIndex",
                       "cannot create index for unopened reader");
        return false;
    }

    // create new index of requested type
    BamIndex* newIndex = BamIndexFactory::CreateIndexOfType(type, reader);
    if ( newIndex == 0 ) {
        stringstream s("");
        s << "could not create index of type: " << type;
        SetErrorString("BamRandomAccessController::CreateIndex", s.str());
        return false;
    }

    // attempt to build index from current BamReader file
    if ( !newIndex->Create() ) {
        const string indexError = newIndex->GetErrorString();
        const string message = "could not create index: \n\t" + indexError;
        SetErrorString("BamRandomAccessController::CreateIndex", message);
        return false;
    }

    // save new index & return success
    SetIndex(newIndex);
    return true;
}

string BamRandomAccessController::GetErrorString(void) const {
    return m_errorString;
}

bool BamRandomAccessController::HasIndex(void) const {
    return ( m_index != 0 );
}

bool BamRandomAccessController::HasRegion(void) const  {
    return ( !m_region.isNull() );
}

bool BamRandomAccessController::IndexHasAlignmentsForReference(const int& refId) {
    return m_index->HasAlignments(refId);
}

bool BamRandomAccessController::LocateIndex(BamReaderPrivate* reader,
                                            const BamIndex::IndexType& preferredType)
{
    // look up index filename, deferring to preferredType if possible
    assert(reader);
    const string& indexFilename = BamIndexFactory::FindIndexFilename(reader->Filename(), preferredType);

    // if no index file found (of any type)
    if ( indexFilename.empty() ) {
        const string message = string("could not find index file for:") + reader->Filename();
        SetErrorString("BamRandomAccessController::LocateIndex", message);
        return false;
    }

    // otherwise open & use index file that was found
    return OpenIndex(indexFilename, reader);
}

bool BamRandomAccessController::OpenIndex(const string& indexFilename, BamReaderPrivate* reader) {

    // attempt create new index of type based on filename
    BamIndex* index = BamIndexFactory::CreateIndexFromFilename(indexFilename, reader);
    if ( index == 0 ) {
        const string message = string("could not open index file: ") + indexFilename;
        SetErrorString("BamRandomAccessController::OpenIndex", message);
        return false;
    }

    // attempt to load data from index file
    if ( !index->Load(indexFilename) ) {
        const string indexError = index->GetErrorString();
        const string message = string("could not load index data from file: ") + indexFilename +
                               "\n\t" + indexError;
        SetErrorString("BamRandomAccessController::OpenIndex", message);
        return false;
    }

    // save new index & return success
    SetIndex(index);
    return true;
}

bool BamRandomAccessController::RegionHasAlignments(void) const {
    return m_hasAlignmentsInRegion;
}

void BamRandomAccessController::SetErrorString(const string& where, const string& what) {
    m_errorString = where + ": " + what;
}

void BamRandomAccessController::SetIndex(BamIndex* index) {
    if ( m_index )
        ClearIndex();
    m_index = index;
}

bool BamRandomAccessController::SetRegion(const BamRegion& region, const int& referenceCount) {

    // store region
    m_region = region;

    // cannot jump when no index is available
    if ( !HasIndex() ) {
        SetErrorString("BamRandomAccessController", "cannot jump if no index data available");
        return false;
    }

    // adjust region as necessary to reflect where data actually begins
    AdjustRegion(referenceCount);

    // if no data present, return true
    //   * Not an error, but future attempts to access alignments in this region will not return data
    //     Returning true is useful in a BamMultiReader setting where some BAM files may
    //     lack alignments in regions where other files still have data available.
    if ( !m_hasAlignmentsInRegion )
        return true;

    // return success/failure of jump to specified region,
    //
    //  * Index::Jump() is allowed to modify the m_hasAlignmentsInRegion flag
    //    This covers 'corner case' where a region is requested that lies beyond the last
    //    alignment on a reference. If this occurs, any subsequent calls to GetNextAlignment[Core]
    //    will not return data. BamMultiReader will still be able to successfully pull alignments
    //    from a region from other files even if this one has no data.
    if ( !m_index->Jump(m_region, &m_hasAlignmentsInRegion) ) {
        const string indexError = m_index->GetErrorString();
        const string message = string("could not set region\n\t") + indexError;
        SetErrorString("BamRandomAccessController::OpenIndex", message);
        return false;
    }
    else
        return true;
}
