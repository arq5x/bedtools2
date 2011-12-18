// ***************************************************************************
// SamSequence.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides direct read/write access to the SAM sequence data fields.
// ***************************************************************************

#include "api/SamSequence.h"
#include <sstream>
using namespace BamTools;
using namespace std;

/*! \struct BamTools::SamSequence
    \brief Represents a SAM sequence entry.

    Provides direct read/write access to the SAM sequence data fields.

    \sa \samSpecURL
*/
/*! \var SamSequence::AssemblyID
    \brief corresponds to \@SQ AS:\<AssemblyID\>
*/
/*! \var SamSequence::Checksum
    \brief corresponds to \@SQ M5:\<Checksum\>
*/
/*! \var SamSequence::Length
    \brief corresponds to \@SQ LN:\<Length\>

    Required for valid SAM header.
*/
/*! \var SamSequence::Name
    \brief corresponds to \@SQ SN:\<Name\>

    Required for valid SAM header.
*/
/*! \var SamSequence::Species
    \brief corresponds to \@SQ SP:\<Species\>
*/
/*! \var SamSequence::URI
    \brief corresponds to \@SQ UR:\<URI\>
*/

/*! \fn SamSequence::SamSequence(void)
    \brief default constructor
*/
SamSequence::SamSequence(void)
    : AssemblyID("")
    , Checksum("")
    , Length("")
    , Name("")
    , Species("")
    , URI("")
{ }

/*! \fn SamSequence::SamSequence(const std::string& name, const int& length)
    \brief constructs sequence with \a name and \a length

    \param name   desired sequence name
    \param length desired sequence length (numeric value)
*/
SamSequence::SamSequence(const std::string& name,
                         const int& length)
    : AssemblyID("")
    , Checksum("")
    , Name(name)
    , Species("")
    , URI("")
{
    stringstream s("");
    s << length;
    Length = s.str();
}

/*! \fn SamSequence::SamSequence(const std::string& name, const std::string& length)
    \brief constructs sequence with \a name and \a length

    \param name   desired sequence name
    \param length desired sequence length (string value)
*/
SamSequence::SamSequence(const std::string& name,
                         const std::string& length)
    : AssemblyID("")
    , Checksum("")
    , Length(length)
    , Name(name)
    , Species("")
    , URI("")
{ }

/*! \fn SamSequence::SamSequence(const SamSequence& other)
    \brief copy constructor
*/
SamSequence::SamSequence(const SamSequence& other)
    : AssemblyID(other.AssemblyID)
    , Checksum(other.Checksum)
    , Length(other.Length)
    , Name(other.Name)
    , Species(other.Species)
    , URI(other.URI)
{ }

/*! \fn SamSequence::~SamSequence(void)
    \brief destructor
*/
SamSequence::~SamSequence(void) { }

/*! \fn void SamSequence::Clear(void)
    \brief Clears all data fields.
*/
void SamSequence::Clear(void) {
    AssemblyID.clear();
    Checksum.clear();
    Length.clear();
    Name.clear();
    Species.clear();
    URI.clear();
}

/*! \fn bool SamSequence::HasAssemblyID(void) const
    \brief Returns \c true if sequence contains \@SQ AS:\<AssemblyID\>
*/
bool SamSequence::HasAssemblyID(void) const {
    return (!AssemblyID.empty());
}

/*! \fn bool SamSequence::HasChecksum(void) const
    \brief Returns \c true if sequence contains \@SQ M5:\<Checksum\>
*/
bool SamSequence::HasChecksum(void) const {
    return (!Checksum.empty());
}

/*! \fn bool SamSequence::HasLength(void) const
    \brief Returns \c true if sequence contains \@SQ LN:\<Length\>
*/
bool SamSequence::HasLength(void) const {
    return (!Length.empty());
}

/*! \fn bool SamSequence::HasName(void) const
    \brief Returns \c true if sequence contains \@SQ SN:\<Name\>
*/
bool SamSequence::HasName(void) const {
    return (!Name.empty());
}

/*! \fn bool SamSequence::HasSpecies(void) const
    \brief Returns \c true if sequence contains \@SQ SP:\<Species\>
*/
bool SamSequence::HasSpecies(void) const {
    return (!Species.empty());
}

/*! \fn bool SamSequence::HasURI(void) const
    \brief Returns \c true if sequence contains \@SQ UR:\<URI\>
*/
bool SamSequence::HasURI(void) const {
    return (!URI.empty());
}
