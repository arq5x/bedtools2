// ***************************************************************************
// SamReadGroup.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides direct read/write access to the SAM read group data fields.
// ***************************************************************************

#include "api/SamReadGroup.h"
using namespace BamTools;
using namespace std;

/*! \struct BamTools::SamReadGroup
    \brief Represents a SAM read group entry.

    Provides direct read/write access to the SAM read group data fields.

    \sa \samSpecURL
*/
/*! \var SamReadGroup::Description
    \brief corresponds to \@RG DS:\<Description\>
*/
/*! \var SamReadGroup::FlowOrder
    \brief corresponds to \@RG FO:\<FlowOrder\>
*/
/*! \var SamReadGroup::ID
    \brief corresponds to \@RG ID:\<ID\>

    Required for valid SAM header.
*/
/*! \var SamReadGroup::KeySequence
    \brief corresponds to \@RG KS:\<KeySequence\>
*/
/*! \var SamReadGroup::Library
    \brief corresponds to \@RG LB:\<Library\>
*/
/*! \var SamReadGroup::PlatformUnit
    \brief corresponds to \@RG PU:\<PlatformUnit\>
*/
/*! \var SamReadGroup::PredictedInsertSize
    \brief corresponds to \@RG PI:\<PredictedInsertSize\>
*/
/*! \var SamReadGroup::ProductionDate
    \brief corresponds to \@RG DT:\<ProductionDate\>
*/
/*! \var SamReadGroup::Program
    \brief corresponds to \@RG PG:\<Program\>
*/
/*! \var SamReadGroup::Sample
    \brief corresponds to \@RG SM:\<Sample\>
*/
/*! \var SamReadGroup::SequencingCenter
    \brief corresponds to \@RG CN:\<SequencingCenter\>
*/
/*! \var SamReadGroup::SequencingTechnology
    \brief corresponds to \@RG PL:\<SequencingTechnology\>
*/

/*! \fn SamReadGroup::SamReadGroup(void)
    \brief default constructor
*/
SamReadGroup::SamReadGroup(void)
    : Description("")
    , FlowOrder("")
    , ID("")
    , KeySequence("")
    , Library("")
    , PlatformUnit("")
    , PredictedInsertSize("")
    , ProductionDate("")
    , Program("")
    , Sample("")
    , SequencingCenter("")
    , SequencingTechnology("")
{ }

/*! \fn SamReadGroup::SamReadGroup(const std::string& id)
    \brief constructs read group with \a id

    \param id desired read group ID
*/
SamReadGroup::SamReadGroup(const std::string& id)
    : Description("")
    , FlowOrder("")
    , ID(id)
    , KeySequence("")
    , Library("")
    , PlatformUnit("")
    , PredictedInsertSize("")
    , ProductionDate("")
    , Program("")
    , Sample("")
    , SequencingCenter("")
    , SequencingTechnology("")
{ }

/*! \fn SamReadGroup::SamReadGroup(const SamReadGroup& other)
    \brief copy constructor
*/
SamReadGroup::SamReadGroup(const SamReadGroup& other)
    : Description(other.Description)
    , FlowOrder(other.FlowOrder)
    , ID(other.ID)
    , KeySequence(other.KeySequence)
    , Library(other.Library)
    , PlatformUnit(other.PlatformUnit)
    , PredictedInsertSize(other.PredictedInsertSize)
    , ProductionDate(other.ProductionDate)
    , Program(other.Program)
    , Sample(other.Sample)
    , SequencingCenter(other.SequencingCenter)
    , SequencingTechnology(other.SequencingTechnology)
{ }

/*! \fn SamReadGroup::~SamReadGroup(void)
    \brief destructor
*/
SamReadGroup::~SamReadGroup(void) { }

/*! \fn void SamReadGroup::Clear(void)
    \brief Clears all data fields.
*/
void SamReadGroup::Clear(void) {
    Description.clear();
    FlowOrder.clear();
    ID.clear();
    KeySequence.clear();
    Library.clear();
    PlatformUnit.clear();
    PredictedInsertSize.clear();
    ProductionDate.clear();
    Program.clear();
    Sample.clear();
    SequencingCenter.clear();
    SequencingTechnology.clear();
}

/*! \fn bool SamReadGroup::HasDescription(void) const
    \brief Returns \c true if read group contains \@RG DS:\<Description\>
*/
bool SamReadGroup::HasDescription(void) const {
    return (!Description.empty());
}

/*! \fn bool SamReadGroup::HasFlowOrder(void) const
    \brief Returns \c true if read group contains \@RG FO:\<FlowOrder\>
*/
bool SamReadGroup::HasFlowOrder(void) const {
    return (!FlowOrder.empty());
}

/*! \fn bool SamReadGroup::HasID(void) const
    \brief Returns \c true if read group contains \@RG: ID:\<ID\>
*/
bool SamReadGroup::HasID(void) const {
    return (!ID.empty());
}

/*! \fn bool SamReadGroup::HasKeySequence(void) const
    \brief Returns \c true if read group contains \@RG KS:\<KeySequence\>
*/
bool SamReadGroup::HasKeySequence(void) const {
    return (!KeySequence.empty());
}

/*! \fn bool SamReadGroup::HasLibrary(void) const
    \brief Returns \c true if read group contains \@RG LB:\<Library\>
*/
bool SamReadGroup::HasLibrary(void) const {
    return (!Library.empty());
}

/*! \fn bool SamReadGroup::HasPlatformUnit(void) const
    \brief Returns \c true if read group contains \@RG PU:\<PlatformUnit\>
*/
bool SamReadGroup::HasPlatformUnit(void) const {
    return (!PlatformUnit.empty());
}

/*! \fn bool SamReadGroup::HasPredictedInsertSize(void) const
    \brief Returns \c true if read group contains \@RG PI:\<PredictedInsertSize\>
*/
bool SamReadGroup::HasPredictedInsertSize(void) const {
    return (!PredictedInsertSize.empty());
}

/*! \fn bool SamReadGroup::HasProductionDate(void) const
    \brief Returns \c true if read group contains \@RG DT:\<ProductionDate\>
*/
bool SamReadGroup::HasProductionDate(void) const {
    return (!ProductionDate.empty());
}

/*! \fn bool SamReadGroup::HasProgram(void) const
    \brief Returns \c true if read group contains \@RG PG:\<Program\>
*/
bool SamReadGroup::HasProgram(void) const {
    return (!Program.empty());
}

/*! \fn bool SamReadGroup::HasSample(void) const
    \brief Returns \c true if read group contains \@RG SM:\<Sample\>
*/
bool SamReadGroup::HasSample(void) const {
    return (!Sample.empty());
}

/*! \fn bool SamReadGroup::HasSequencingCenter(void) const
    \brief Returns \c true if read group contains \@RG CN:\<SequencingCenter\>
*/
bool SamReadGroup::HasSequencingCenter(void) const {
    return (!SequencingCenter.empty());
}

/*! \fn bool SamReadGroup::HasSequencingTechnology(void) const
    \brief Returns \c true if read group contains \@RG PL:\<SequencingTechnology\>
*/
bool SamReadGroup::HasSequencingTechnology(void) const {
    return (!SequencingTechnology.empty());
}
