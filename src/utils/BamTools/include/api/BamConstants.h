// ***************************************************************************
// BamConstants.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides basic constants for handling BAM files.
// ***************************************************************************

#ifndef BAM_CONSTANTS_H
#define BAM_CONSTANTS_H

#include <string>

/*! \namespace BamTools::Constants
    \brief Provides basic constants for handling BAM files.
*/

namespace BamTools {
namespace Constants {

const int BAM_SIZEOF_INT = 4;

// header magic number
const char* const  BAM_HEADER_MAGIC = "BAM\1";
const unsigned int BAM_HEADER_MAGIC_LENGTH = 4;

// BAM alignment core size
const int BAM_CORE_SIZE = 32;
const int BAM_CORE_BUFFER_SIZE = 8;

// BAM alignment flags
const int BAM_ALIGNMENT_PAIRED              = 0x0001;
const int BAM_ALIGNMENT_PROPER_PAIR         = 0x0002;
const int BAM_ALIGNMENT_UNMAPPED            = 0x0004;
const int BAM_ALIGNMENT_MATE_UNMAPPED       = 0x0008;
const int BAM_ALIGNMENT_REVERSE_STRAND      = 0x0010;
const int BAM_ALIGNMENT_MATE_REVERSE_STRAND = 0x0020;
const int BAM_ALIGNMENT_READ_1              = 0x0040;
const int BAM_ALIGNMENT_READ_2              = 0x0080;
const int BAM_ALIGNMENT_SECONDARY           = 0x0100;
const int BAM_ALIGNMENT_QC_FAILED           = 0x0200;
const int BAM_ALIGNMENT_DUPLICATE           = 0x0400;

// CIGAR constants
const char* const BAM_CIGAR_LOOKUP = "MIDNSHP=X";
const int BAM_CIGAR_MATCH    = 0;
const int BAM_CIGAR_INS      = 1;
const int BAM_CIGAR_DEL      = 2;
const int BAM_CIGAR_REFSKIP  = 3;
const int BAM_CIGAR_SOFTCLIP = 4;
const int BAM_CIGAR_HARDCLIP = 5;
const int BAM_CIGAR_PAD      = 6;
const int BAM_CIGAR_SEQMATCH = 7;
const int BAM_CIGAR_MISMATCH = 8;

const char BAM_CIGAR_MATCH_CHAR    = 'M';
const char BAM_CIGAR_INS_CHAR      = 'I';
const char BAM_CIGAR_DEL_CHAR      = 'D';
const char BAM_CIGAR_REFSKIP_CHAR  = 'N';
const char BAM_CIGAR_SOFTCLIP_CHAR = 'S';
const char BAM_CIGAR_HARDCLIP_CHAR = 'H';
const char BAM_CIGAR_PAD_CHAR      = 'P';
const char BAM_CIGAR_SEQMATCH_CHAR = '=';
const char BAM_CIGAR_MISMATCH_CHAR = 'X';

const int BAM_CIGAR_SHIFT    = 4;
const int BAM_CIGAR_MASK     = ((1 << BAM_CIGAR_SHIFT) - 1);

// BAM tag types
const char BAM_TAG_TYPE_ASCII  = 'A';
const char BAM_TAG_TYPE_UINT8  = 'c';
const char BAM_TAG_TYPE_INT8   = 'C';
const char BAM_TAG_TYPE_UINT16 = 's';
const char BAM_TAG_TYPE_INT16  = 'S';
const char BAM_TAG_TYPE_UINT32 = 'i';
const char BAM_TAG_TYPE_INT32  = 'I';
const char BAM_TAG_TYPE_FLOAT  = 'f';
const char BAM_TAG_TYPE_STRING = 'Z';
const char BAM_TAG_TYPE_HEX    = 'H';
const char BAM_TAG_TYPE_ARRAY  = 'B';

const size_t BAM_TAG_TAGSIZE  = 2;
const size_t BAM_TAG_TYPESIZE = 1;
const int BAM_TAG_ARRAYBASE_SIZE = 8;

// DNA bases
const char* const BAM_DNA_LOOKUP = "=ACMGRSVTWYHKDBN";
const unsigned char BAM_BASECODE_EQUAL = 0;
const unsigned char BAM_BASECODE_A     = 1;
const unsigned char BAM_BASECODE_C     = 2;
const unsigned char BAM_BASECODE_G     = 4;
const unsigned char BAM_BASECODE_T     = 8;
const unsigned char BAM_BASECODE_N     = 15;

const char BAM_DNA_EQUAL   = '=';
const char BAM_DNA_A       = 'A';
const char BAM_DNA_C       = 'C';
const char BAM_DNA_G       = 'G';
const char BAM_DNA_T       = 'T';
const char BAM_DNA_N       = 'N';
const char BAM_DNA_DEL     = '-';
const char BAM_DNA_PAD     = '*';

// zlib constants
const int GZIP_ID1   = 31;
const int GZIP_ID2   = 139;
const int CM_DEFLATE = 8;
const int FLG_FEXTRA = 4;
const int OS_UNKNOWN = 255;
const int BGZF_XLEN  = 6;
const int BGZF_ID1   = 66;
const int BGZF_ID2   = 67;
const int BGZF_LEN   = 2;
const int GZIP_WINDOW_BITS    = -15;
const int Z_DEFAULT_MEM_LEVEL = 8;

// BZGF constants
const int BGZF_BLOCK_HEADER_LENGTH = 18;
const int BGZF_BLOCK_FOOTER_LENGTH = 8;
const int BGZF_MAX_BLOCK_SIZE      = 65536;
const int BGZF_DEFAULT_BLOCK_SIZE  = 65536;

} // namespace Constants
} // namespace BamTools

#endif // BAM_CONSTANTS_H
