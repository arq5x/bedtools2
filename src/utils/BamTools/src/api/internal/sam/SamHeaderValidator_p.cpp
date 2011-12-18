// ***************************************************************************
// SamHeaderValidator.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides functionality for validating SamHeader data
// ***************************************************************************

#include "api/SamConstants.h"
#include "api/SamHeader.h"
#include "api/internal/sam/SamHeaderValidator_p.h"
#include "api/internal/sam/SamHeaderVersion_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cctype>
#include <set>
#include <sstream>
using namespace std;

// ------------------------
// static utility methods
// -------------------------

static
bool caseInsensitiveCompare(const string& lhs, const string& rhs) {

    // can omit checking chars if lengths not equal
    const int lhsLength = lhs.length();
    const int rhsLength = rhs.length();
    if ( lhsLength != rhsLength )
        return false;

    // do *basic* toupper checks on each string char's
    for ( int i = 0; i < lhsLength; ++i ) {
        if ( toupper( (int)lhs.at(i)) != toupper( (int)rhs.at(i)) )
            return false;
    }

    // otherwise OK
    return true;
}

// ------------------------------------------------------------------------
// Allow validation rules to vary, as needed, between SAM header versions
//
// use SAM_VERSION_X_Y to tag important changes
//
// Together, they will allow for comparisons like:
// if ( m_version < SAM_VERSION_2_0 ) {
//     // use some older rule
// else
//     // use rule introduced with version 2.0

static const SamHeaderVersion SAM_VERSION_1_0 = SamHeaderVersion(1,0);
static const SamHeaderVersion SAM_VERSION_1_1 = SamHeaderVersion(1,1);
static const SamHeaderVersion SAM_VERSION_1_2 = SamHeaderVersion(1,2);
static const SamHeaderVersion SAM_VERSION_1_3 = SamHeaderVersion(1,3);
static const SamHeaderVersion SAM_VERSION_1_4 = SamHeaderVersion(1,4);

// TODO: This functionality is currently unused.
//       Make validation "version-aware."
//
// ------------------------------------------------------------------------

const string SamHeaderValidator::ERROR_PREFIX = "ERROR: ";
const string SamHeaderValidator::WARN_PREFIX  = "WARNING: ";
const string SamHeaderValidator::NEWLINE      = "\n";

SamHeaderValidator::SamHeaderValidator(const SamHeader& header)
    : m_header(header)
{ }

SamHeaderValidator::~SamHeaderValidator(void) { }

void SamHeaderValidator::AddError(const string& message) {
    m_errorMessages.push_back(ERROR_PREFIX + message + NEWLINE);
}

void SamHeaderValidator::AddWarning(const string& message) {
    m_warningMessages.push_back(WARN_PREFIX + message + NEWLINE);
}

void SamHeaderValidator::PrintErrorMessages(ostream& stream) {

    // skip if no error messages
    if ( m_errorMessages.empty() )
        return;

    // print error header line
    stream << "* SAM header has " << m_errorMessages.size() << " errors:" << endl;

    // print each error message
    vector<string>::const_iterator errorIter = m_errorMessages.begin();
    vector<string>::const_iterator errorEnd  = m_errorMessages.end();
    for ( ; errorIter != errorEnd; ++errorIter )
        stream << (*errorIter);
}

void SamHeaderValidator::PrintMessages(ostream& stream) {
    PrintErrorMessages(stream);
    PrintWarningMessages(stream);
}

void SamHeaderValidator::PrintWarningMessages(ostream& stream) {

    // skip if no warning messages
    if ( m_warningMessages.empty() )
        return;

    // print warning header line
    stream << "* SAM header has " << m_warningMessages.size() << " warnings:" << endl;

    // print each warning message
    vector<string>::const_iterator warnIter = m_warningMessages.begin();
    vector<string>::const_iterator warnEnd  = m_warningMessages.end();
    for ( ; warnIter != warnEnd; ++warnIter )
        stream << (*warnIter);
}

// entry point for validation
bool SamHeaderValidator::Validate(void) {
    bool isValid = true;
    isValid &= ValidateMetadata();
    isValid &= ValidateSequenceDictionary();
    isValid &= ValidateReadGroupDictionary();
    isValid &= ValidateProgramChain();
    return isValid;
}

// check all SAM header 'metadata'
bool SamHeaderValidator::ValidateMetadata(void) {
    bool isValid = true;
    isValid &= ValidateVersion();
    isValid &= ValidateSortOrder();
    isValid &= ValidateGroupOrder();
    return isValid;
}

// check SAM header version tag
bool SamHeaderValidator::ValidateVersion(void) {

    const string& version = m_header.Version;

    // warn if version not present
    if ( version.empty() ) {
        AddWarning("Version (VN) missing. Not required, but strongly recommended");
        return true;
    }

    // invalid if version does not contain a period
    const size_t periodFound = version.find(Constants::SAM_PERIOD);
    if ( periodFound == string::npos ) {
        AddError("Invalid version (VN) format: " + version);
        return false;
    }

    // invalid if major version is empty or contains non-digits
    const string majorVersion = version.substr(0, periodFound);
    if ( majorVersion.empty() || !ContainsOnlyDigits(majorVersion) ) {
        AddError("Invalid version (VN) format: " + version);
        return false;
    }

    // invalid if major version is empty or contains non-digits
    const string minorVersion = version.substr(periodFound + 1);
    if ( minorVersion.empty() || !ContainsOnlyDigits(minorVersion) ) {
        AddError("Invalid version (VN) format: " + version);
        return false;
    }

    // TODO: check if version is not just syntactically OK,
    // but is also a valid SAM version ( 1.0 .. CURRENT )

    // all checked out this far, then version is OK
    return true;
}

// assumes non-empty input string
bool SamHeaderValidator::ContainsOnlyDigits(const string& s) {
    const size_t nonDigitPosition = s.find_first_not_of(Constants::SAM_DIGITS);
    return ( nonDigitPosition == string::npos ) ;
}

// validate SAM header sort order tag
bool SamHeaderValidator::ValidateSortOrder(void) {

    const string& sortOrder = m_header.SortOrder;

    // warn if sort order not present
    if ( sortOrder.empty() ) {
        AddWarning("Sort order (SO) missing. Not required, but strongly recommended");
        return true;
    }

    // if sort order is valid keyword
    if ( sortOrder == Constants::SAM_HD_SORTORDER_COORDINATE ||
         sortOrder == Constants::SAM_HD_SORTORDER_QUERYNAME  ||
         sortOrder == Constants::SAM_HD_SORTORDER_UNSORTED
       )
    {
        return true;
    }

    // otherwise
    AddError("Invalid sort order (SO): " + sortOrder);
    return false;
}

// validate SAM header group order tag
bool SamHeaderValidator::ValidateGroupOrder(void) {

    const string& groupOrder = m_header.GroupOrder;

    // if no group order, no problem, just return OK
    if ( groupOrder.empty() )
        return true;

    // if group order is valid keyword
    if ( groupOrder == Constants::SAM_HD_GROUPORDER_NONE  ||
         groupOrder == Constants::SAM_HD_GROUPORDER_QUERY ||
         groupOrder == Constants::SAM_HD_GROUPORDER_REFERENCE
       )
    {
        return true;
    }

    // otherwise
    AddError("Invalid group order (GO): " + groupOrder);
    return false;
}

// validate SAM header sequence dictionary
bool SamHeaderValidator::ValidateSequenceDictionary(void) {

    bool isValid = true;

    // check for unique sequence names
    isValid &= ContainsUniqueSequenceNames();

    // iterate over sequences
    const SamSequenceDictionary& sequences = m_header.Sequences;
    SamSequenceConstIterator seqIter = sequences.ConstBegin();
    SamSequenceConstIterator seqEnd  = sequences.ConstEnd();
    for ( ; seqIter != seqEnd; ++seqIter ) {
        const SamSequence& seq = (*seqIter);
        isValid &= ValidateSequence(seq);
    }

    // return validation state
    return isValid;
}

// make sure all SQ names are unique
bool SamHeaderValidator::ContainsUniqueSequenceNames(void) {

    bool isValid = true;
    set<string> sequenceNames;
    set<string>::iterator nameIter;

    // iterate over sequences
    const SamSequenceDictionary& sequences = m_header.Sequences;
    SamSequenceConstIterator seqIter = sequences.ConstBegin();
    SamSequenceConstIterator seqEnd  = sequences.ConstEnd();
    for ( ; seqIter != seqEnd; ++seqIter ) {
        const SamSequence& seq = (*seqIter);

        // lookup sequence name
        const string& name = seq.Name;
        nameIter = sequenceNames.find(name);

        // error if found (duplicate entry)
        if ( nameIter != sequenceNames.end() ) {
            AddError("Sequence name (SN): " + name + " is not unique");
            isValid = false;
        }

        // otherwise ok, store name
        sequenceNames.insert(name);
    }

    // return validation state
    return isValid;
}

// validate SAM header sequence entry
bool SamHeaderValidator::ValidateSequence(const SamSequence& seq) {
    bool isValid = true;
    isValid &= CheckNameFormat(seq.Name);
    isValid &= CheckLengthInRange(seq.Length);
    return isValid;
}

// check sequence name is valid format
bool SamHeaderValidator::CheckNameFormat(const string& name) {

    // invalid if name is empty
    if ( name.empty() ) {
        AddError("Sequence entry (@SQ) is missing SN tag");
        return false;
    }

    // invalid if first character is a reserved char
    const char firstChar = name.at(0);
    if ( firstChar == Constants::SAM_EQUAL || firstChar == Constants::SAM_STAR ) {
        AddError("Invalid sequence name (SN): " + name);
        return false;
    }
    // otherwise OK
    return true;
}

// check that sequence length is within accepted range
bool SamHeaderValidator::CheckLengthInRange(const string& length) {

    // invalid if empty
    if ( length.empty() ) {
        AddError("Sequence entry (@SQ) is missing LN tag");
        return false;
    }

    // convert string length to numeric
    stringstream lengthStream(length);
    unsigned int sequenceLength;
    lengthStream >> sequenceLength;

    // invalid if length outside accepted range
    if ( sequenceLength < Constants::SAM_SQ_LENGTH_MIN || sequenceLength > Constants::SAM_SQ_LENGTH_MAX ) {
        AddError("Sequence length (LN): " + length + " out of range");
        return false;
    }

    // otherwise OK
    return true;
}

// validate SAM header read group dictionary
bool SamHeaderValidator::ValidateReadGroupDictionary(void) {

    bool isValid = true;

    // check for unique read group IDs & platform units
    isValid &= ContainsUniqueIDsAndPlatformUnits();

    // iterate over read groups
    const SamReadGroupDictionary& readGroups = m_header.ReadGroups;
    SamReadGroupConstIterator rgIter = readGroups.ConstBegin();
    SamReadGroupConstIterator rgEnd  = readGroups.ConstEnd();
    for ( ; rgIter != rgEnd; ++rgIter ) {
        const SamReadGroup& rg = (*rgIter);
        isValid &= ValidateReadGroup(rg);
    }

    // return validation state
    return isValid;
}

// make sure RG IDs and platform units are unique
bool SamHeaderValidator::ContainsUniqueIDsAndPlatformUnits(void) {

    bool isValid = true;
    set<string> readGroupIds;
    set<string> platformUnits;
    set<string>::iterator idIter;
    set<string>::iterator puIter;

    // iterate over sequences
    const SamReadGroupDictionary& readGroups = m_header.ReadGroups;
    SamReadGroupConstIterator rgIter = readGroups.ConstBegin();
    SamReadGroupConstIterator rgEnd  = readGroups.ConstEnd();
    for ( ; rgIter != rgEnd; ++rgIter ) {
        const SamReadGroup& rg = (*rgIter);

        // --------------------------------
        // check for unique ID

        // lookup read group ID
        const string& id = rg.ID;
        idIter = readGroupIds.find(id);

        // error if found (duplicate entry)
        if ( idIter != readGroupIds.end() ) {
            AddError("Read group ID (ID): " + id + " is not unique");
            isValid = false;
        }

        // otherwise ok, store id
        readGroupIds.insert(id);

        // --------------------------------
        // check for unique platform unit

        // lookup platform unit
        const string& pu = rg.PlatformUnit;
        puIter = platformUnits.find(pu);

        // error if found (duplicate entry)
        if ( puIter != platformUnits.end() ) {
            AddError("Platform unit (PU): " + pu + " is not unique");
            isValid = false;
        }

        // otherwise ok, store platform unit
        platformUnits.insert(pu);
    }

    // return validation state
    return isValid;
}

// validate SAM header read group entry
bool SamHeaderValidator::ValidateReadGroup(const SamReadGroup& rg) {
    bool isValid = true;
    isValid &= CheckReadGroupID(rg.ID);
    isValid &= CheckSequencingTechnology(rg.SequencingTechnology);
    return isValid;
}

// make sure RG ID exists
bool SamHeaderValidator::CheckReadGroupID(const string& id) {

    // invalid if empty
    if ( id.empty() ) {
        AddError("Read group entry (@RG) is missing ID tag");
        return false;
    }

    // otherwise OK
    return true;
}

// make sure RG sequencing tech is one of the accepted keywords
bool SamHeaderValidator::CheckSequencingTechnology(const string& technology) {

    // if no technology provided, no problem, just return OK
    if ( technology.empty() )
        return true;

    // if technology is valid keyword
    if ( caseInsensitiveCompare(technology, Constants::SAM_RG_SEQTECHNOLOGY_CAPILLARY)  ||
         caseInsensitiveCompare(technology, Constants::SAM_RG_SEQTECHNOLOGY_HELICOS)    ||
         caseInsensitiveCompare(technology, Constants::SAM_RG_SEQTECHNOLOGY_ILLUMINA)   ||
         caseInsensitiveCompare(technology, Constants::SAM_RG_SEQTECHNOLOGY_IONTORRENT) ||
         caseInsensitiveCompare(technology, Constants::SAM_RG_SEQTECHNOLOGY_LS454)      ||
         caseInsensitiveCompare(technology, Constants::SAM_RG_SEQTECHNOLOGY_PACBIO)     ||
         caseInsensitiveCompare(technology, Constants::SAM_RG_SEQTECHNOLOGY_SOLID)
       )
    {
        return true;
    }

    // otherwise
    AddError("Invalid read group sequencing platform (PL): " + technology);
    return false;
}

// validate the SAM header "program chain"
bool SamHeaderValidator::ValidateProgramChain(void) {
    bool isValid = true;
    isValid &= ContainsUniqueProgramIds();
    isValid &= ValidatePreviousProgramIds();
    return isValid;
}

// make sure all PG IDs are unique
bool SamHeaderValidator::ContainsUniqueProgramIds(void) {

    bool isValid = true;
    set<string> programIds;
    set<string>::iterator pgIdIter;

    // iterate over program records
    const SamProgramChain& programs = m_header.Programs;
    SamProgramConstIterator pgIter = programs.ConstBegin();
    SamProgramConstIterator pgEnd  = programs.ConstEnd();
    for ( ; pgIter != pgEnd; ++pgIter ) {
        const SamProgram& pg = (*pgIter);

        // lookup program ID
        const string& pgId = pg.ID;
        pgIdIter = programIds.find(pgId);

        // error if found (duplicate entry)
        if ( pgIdIter != programIds.end() ) {
            AddError("Program ID (ID): " + pgId + " is not unique");
            isValid = false;
        }

        // otherwise ok, store ID
        programIds.insert(pgId);
    }

    // return validation state
    return isValid;
}

// make sure that any PP tags present point to existing @PG IDs
bool SamHeaderValidator::ValidatePreviousProgramIds(void) {

    bool isValid = true;

    // iterate over program records
    const SamProgramChain& programs = m_header.Programs;
    SamProgramConstIterator pgIter = programs.ConstBegin();
    SamProgramConstIterator pgEnd  = programs.ConstEnd();
    for ( ; pgIter != pgEnd; ++pgIter ) {
        const SamProgram& pg = (*pgIter);

        // ignore record for validation if PreviousProgramID is empty
        const string& ppId = pg.PreviousProgramID;
        if ( ppId.empty() )
            continue;

        // see if program "chain" contains an entry for ppId
        if ( !programs.Contains(ppId) ) {
            AddError("PreviousProgramID (PP): " + ppId + " is not a known ID");
            isValid = false;
        }
    }

    // return validation state
    return isValid;
}
