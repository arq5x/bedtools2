// ***************************************************************************
// SamHeader.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides direct read/write access to the SAM header data fields.
// ***************************************************************************

#include "api/SamConstants.h"
#include "api/SamHeader.h"
#include "api/internal/utils/BamException_p.h"
#include "api/internal/sam/SamFormatParser_p.h"
#include "api/internal/sam/SamFormatPrinter_p.h"
#include "api/internal/sam/SamHeaderValidator_p.h"
using namespace BamTools;
using namespace BamTools::Internal;
using namespace std;

/*! \struct BamTools::SamHeader
    \brief Represents the SAM-formatted text header that is part of the BAM file header.

    Provides direct read/write access to the SAM header data fields.

    \sa \samSpecURL
*/
/*! \var SamHeader::Version
    \brief corresponds to \@HD VN:\<Version\>

    Required for valid SAM header, if \@HD record is present.
*/
/*! \var SamHeader::SortOrder
    \brief corresponds to \@HD SO:\<SortOrder\>
*/
/*! \var SamHeader::GroupOrder
    \brief corresponds to \@HD GO:\<GroupOrder\>
*/
/*! \var SamHeader::Sequences
    \brief corresponds to \@SQ entries
    \sa SamSequence, SamSequenceDictionary
*/
/*! \var SamHeader::ReadGroups
    \brief corresponds to \@RG entries
    \sa SamReadGroup, SamReadGroupDictionary
*/
/*! \var SamHeader::Programs
    \brief corresponds to \@PG entries
    \sa SamProgram, SamProgramChain
*/
/*! \var SamHeader::Comments
    \brief corresponds to \@CO entries
*/

/*! \fn SamHeader::SamHeader(const std::string& headerText = "")
    \brief constructor
*/
SamHeader::SamHeader(const std::string& headerText)
    : Version("")
    , SortOrder(Constants::SAM_HD_SORTORDER_UNKNOWN)
    , GroupOrder("")
{
    SetHeaderText(headerText);
}

/*! \fn SamHeader::SamHeader(const SamHeader& other)
    \brief copy constructor
*/
SamHeader::SamHeader(const SamHeader& other)
    : Version(other.Version)
    , SortOrder(other.SortOrder)
    , GroupOrder(other.GroupOrder)
    , Sequences(other.Sequences)
    , ReadGroups(other.ReadGroups)
    , Programs(other.Programs)
{ }

/*! \fn SamHeader::~SamHeader(void)
    \brief destructor
*/
SamHeader::~SamHeader(void) { }

/*! \fn void SamHeader::Clear(void)
    \brief Clears all header contents.
*/
void SamHeader::Clear(void) {

    // clear SAM header components
    Version.clear();
    SortOrder.clear();
    GroupOrder.clear();
    Sequences.Clear();
    ReadGroups.Clear();
    Programs.Clear();
    Comments.clear();

    // clear error string
    m_errorString.clear();
}

/*! \fn std::string SamHeader::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
std::string SamHeader::GetErrorString(void) const {
    return m_errorString;
}

/*! \fn bool SamHeader::HasError(void) const
    \brief Returns \c true if header encountered an error
*/
bool SamHeader::HasError(void) const {
    return (!m_errorString.empty());
}

/*! \fn bool SamHeader::HasVersion(void) const
    \brief Returns \c true if header contains \@HD ID:\<Version\>
*/
bool SamHeader::HasVersion(void) const {
    return (!Version.empty());
}

/*! \fn bool SamHeader::HasSortOrder(void) const
    \brief Returns \c true if header contains \@HD SO:\<SortOrder\>
*/
bool SamHeader::HasSortOrder(void) const {
    return (!SortOrder.empty());
}

/*! \fn bool SamHeader::HasGroupOrder(void) const
    \brief Returns \c true if header contains \@HD GO:\<GroupOrder\>
*/
bool SamHeader::HasGroupOrder(void) const {
    return (!GroupOrder.empty());
}

/*! \fn bool SamHeader::HasSequences(void) const
    \brief Returns \c true if header contains any \@SQ entries
*/
bool SamHeader::HasSequences(void) const {
    return (!Sequences.IsEmpty());
}

/*! \fn bool SamHeader::HasReadGroups(void) const
    \brief Returns \c true if header contains any \@RG entries
*/
bool SamHeader::HasReadGroups(void) const {
    return (!ReadGroups.IsEmpty());
}

/*! \fn bool SamHeader::HasPrograms(void) const
    \brief Returns \c true if header contains any \@PG entries
*/
bool SamHeader::HasPrograms(void) const {
    return (!Programs.IsEmpty());
}

/*! \fn bool SamHeader::HasComments(void) const
    \brief Returns \c true if header contains any \@CO entries
*/
bool SamHeader::HasComments(void) const {
    return (!Comments.empty());
}

/*! \fn bool SamHeader::IsValid(bool verbose = false) const
    \brief Checks header contents for required data and proper formatting.

    \param[in] verbose If set to true, validation errors & warnings will be printed to stderr.
                       Otherwise, messages are available through SamHeader::GetErrorString().
    \return \c true if SAM header is well-formed
*/
bool SamHeader::IsValid(bool verbose) const {

    SamHeaderValidator validator(*this);

    // if SAM header is valid, return success
    if ( validator.Validate() )
        return true;

    // otherwiser
    else {

        // print messages to stderr
        if ( verbose )
            validator.PrintMessages(std::cerr);

        // or catch in local error string
        else {
            stringstream errorStream("");
            validator.PrintMessages(errorStream);
            m_errorString = errorStream.str();
        }
        return false;
    }
}

/*! \fn void SamHeader::SetHeaderText(const std::string& headerText)
    \brief Replaces header contents with \a headerText.

    \param[in] headerText SAM formatted-text that will be parsed into data fields
*/
void SamHeader::SetHeaderText(const std::string& headerText) {

    // clear prior data
    Clear();

    try {
        SamFormatParser parser(*this);
        parser.Parse(headerText);
    } catch ( BamException& e ) {

        // clear anything parsed so far
        // no telling what's valid and what's partially parsed
        Clear();

        // set error string
        m_errorString = e.what();
    }
}

/*! \fn std::string SamHeader::ToString(void) const
    \brief Converts data fields to SAM-formatted text.

    Applies any local modifications made since creating this object or calling SetHeaderText().

    \return SAM-formatted header text
*/
string SamHeader::ToString(void) const {
    SamFormatPrinter printer(*this);
    return printer.ToString();
}
