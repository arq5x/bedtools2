// ***************************************************************************
// SamProgram.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides direct read/write access to the SAM header program records.
// ***************************************************************************

#include "api/SamProgram.h"
using namespace BamTools;
using namespace std;

/*! \struct BamTools::SamProgram
    \brief Represents a SAM program record.

    Provides direct read/write access to the SAM header program records.

    \sa \samSpecURL
*/
/*! \var SamProgram::CommandLine
    \brief corresponds to \@PG CL:\<CommandLine\>
*/
/*! \var SamProgram::ID
    \brief corresponds to \@PG ID:\<ID\>

    Required for valid SAM header.
*/
/*! \var SamProgram::Name
    \brief corresponds to \@PG PN:\<Name\>
*/
/*! \var SamProgram::PreviousProgramID
    \brief corresponds to \@PG PP:\<PreviousProgramID\>
*/
/*! \var SamProgram::Version
    \brief corresponds to \@PG VN:\<Version\>
*/
/*! \var SamProgram::NextProgramID
    \internal
    Holds ID of the "next" program record in a SamProgramChain
*/

/*! \fn SamProgram::SamProgram(void)
    \brief default constructor
*/
SamProgram::SamProgram(void)
    : CommandLine("")
    , ID("")
    , Name("")
    , PreviousProgramID("")
    , Version("")
    , NextProgramID("")
{ }

/*! \fn SamProgram::SamProgram(const std::string& id)
    \brief constructs program record with \a id

    \param id desired program record ID
*/
SamProgram::SamProgram(const std::string& id)
    : CommandLine("")
    , ID(id)
    , Name("")
    , PreviousProgramID("")
    , Version("")
    , NextProgramID("")
{ }

/*! \fn SamProgram::SamProgram(const SamProgram& other)
    \brief copy constructor
*/
SamProgram::SamProgram(const SamProgram& other)
    : CommandLine(other.CommandLine)
    , ID(other.ID)
    , Name(other.Name)
    , PreviousProgramID(other.PreviousProgramID)
    , Version(other.Version)
    , NextProgramID(other.NextProgramID)
{ }

/*! \fn SamProgram::~SamProgram(void)
    \brief destructor
*/
SamProgram::~SamProgram(void) { }

/*! \fn void SamProgram::Clear(void)
    \brief Clears all data fields.
*/
void SamProgram::Clear(void) {
    CommandLine.clear();
    ID.clear();
    Name.clear();
    PreviousProgramID.clear();
    Version.clear();
    NextProgramID.clear();
}

/*! \fn bool SamProgram::HasCommandLine(void) const
    \brief Returns \c true if program record contains \@PG: CL:\<CommandLine\>
*/
bool SamProgram::HasCommandLine(void) const {
    return (!CommandLine.empty());
}

/*! \fn bool SamProgram::HasID(void) const
    \brief Returns \c true if program record contains \@PG: ID:\<ID\>
*/
bool SamProgram::HasID(void) const {
    return (!ID.empty());
}

/*! \fn bool SamProgram::HasName(void) const
    \brief Returns \c true if program record contains \@PG: PN:\<Name\>
*/
bool SamProgram::HasName(void) const {
    return (!Name.empty());
}

/*! \fn bool SamProgram::HasNextProgramID(void) const
    \internal
    \return true if program has a "next" record in a SamProgramChain
*/
bool SamProgram::HasNextProgramID(void) const {
    return (!NextProgramID.empty());
}

/*! \fn bool SamProgram::HasPreviousProgramID(void) const
    \brief Returns \c true if program record contains \@PG: PP:\<PreviousProgramID\>
*/
bool SamProgram::HasPreviousProgramID(void) const {
    return (!PreviousProgramID.empty());
}

/*! \fn bool SamProgram::HasVersion(void) const
    \brief Returns \c true if program record contains \@PG: VN:\<Version\>
*/
bool SamProgram::HasVersion(void) const {
    return (!Version.empty());
}
