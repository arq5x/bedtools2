// ***************************************************************************
// SamConstants.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 27 March 2012 (DB)
// ---------------------------------------------------------------------------
// Provides constants for SAM header
// ***************************************************************************

#ifndef SAM_CONSTANTS_H
#define SAM_CONSTANTS_H

#include <string>

namespace BamTools {
namespace Constants {

// basic char constants used in SAM format
const char SAM_COLON  = ':';
const char SAM_EQUAL  = '=';
const char SAM_PERIOD = '.';
const char SAM_STAR   = '*';
const char SAM_TAB    = '\t';
const std::string SAM_DIGITS = "0123456789";

const std::string SAM_CURRENT_VERSION = "1.4";

// HD entries
const std::string SAM_HD_BEGIN_TOKEN    = "@HD";
const std::string SAM_HD_VERSION_TAG    = "VN";
const std::string SAM_HD_SORTORDER_TAG  = "SO";
const std::string SAM_HD_GROUPORDER_TAG = "GO";

// SQ entries
const std::string SAM_SQ_BEGIN_TOKEN    = "@SQ";
const std::string SAM_SQ_ASSEMBLYID_TAG = "AS";
const std::string SAM_SQ_CHECKSUM_TAG   = "M5";
const std::string SAM_SQ_LENGTH_TAG     = "LN";
const std::string SAM_SQ_NAME_TAG       = "SN";
const std::string SAM_SQ_SPECIES_TAG    = "SP";
const std::string SAM_SQ_URI_TAG        = "UR";

// RG entries
const std::string SAM_RG_BEGIN_TOKEN             = "@RG";
const std::string SAM_RG_DESCRIPTION_TAG         = "DS";
const std::string SAM_RG_FLOWORDER_TAG           = "FO";
const std::string SAM_RG_ID_TAG                  = "ID";
const std::string SAM_RG_KEYSEQUENCE_TAG         = "KS";
const std::string SAM_RG_LIBRARY_TAG             = "LB";
const std::string SAM_RG_PLATFORMUNIT_TAG        = "PU";
const std::string SAM_RG_PREDICTEDINSERTSIZE_TAG = "PI";
const std::string SAM_RG_PRODUCTIONDATE_TAG      = "DT";
const std::string SAM_RG_PROGRAM_TAG             = "PG";
const std::string SAM_RG_SAMPLE_TAG              = "SM";
const std::string SAM_RG_SEQCENTER_TAG           = "CN";
const std::string SAM_RG_SEQTECHNOLOGY_TAG       = "PL";

// PG entries
const std::string SAM_PG_BEGIN_TOKEN         = "@PG";
const std::string SAM_PG_COMMANDLINE_TAG     = "CL";
const std::string SAM_PG_ID_TAG              = "ID";
const std::string SAM_PG_NAME_TAG            = "PN";
const std::string SAM_PG_PREVIOUSPROGRAM_TAG = "PP";
const std::string SAM_PG_VERSION_TAG         = "VN";

// CO entries
const std::string SAM_CO_BEGIN_TOKEN = "@CO";

// HD:SO values
const std::string SAM_HD_SORTORDER_COORDINATE = "coordinate";
const std::string SAM_HD_SORTORDER_QUERYNAME  = "queryname";
const std::string SAM_HD_SORTORDER_UNKNOWN    = "unknown";
const std::string SAM_HD_SORTORDER_UNSORTED   = "unsorted";

// HD:GO values
const std::string SAM_HD_GROUPORDER_NONE      = "none";
const std::string SAM_HD_GROUPORDER_QUERY     = "query";
const std::string SAM_HD_GROUPORDER_REFERENCE = "reference";

// SQ:LN values
const unsigned int SAM_SQ_LENGTH_MIN = 1;
const unsigned int SAM_SQ_LENGTH_MAX = 536870911; // 2^29 - 1

// RG:PL values
const std::string SAM_RG_SEQTECHNOLOGY_CAPILLARY  = "CAPILLARY";
const std::string SAM_RG_SEQTECHNOLOGY_HELICOS    = "HELICOS";
const std::string SAM_RG_SEQTECHNOLOGY_ILLUMINA   = "ILLUMINA";
const std::string SAM_RG_SEQTECHNOLOGY_IONTORRENT = "IONTORRENT";
const std::string SAM_RG_SEQTECHNOLOGY_LS454      = "LS454";
const std::string SAM_RG_SEQTECHNOLOGY_PACBIO     = "PACBIO";
const std::string SAM_RG_SEQTECHNOLOGY_SOLID      = "SOLID";

} // namespace Constants
} // namespace BamTools

#endif // SAM_CONSTANTS_H
