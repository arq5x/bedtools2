// ***************************************************************************
// BamConstants.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 16 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides basic constants for handling BAM files.
// ***************************************************************************

#ifndef BAM_CONSTANTS_H
#define BAM_CONSTANTS_H

#include "api/api_global.h"
#include <cassert>
#include <string>

/*! \namespace BamTools::Constants
    \brief Provides basic constants for handling BAM files.
*/

namespace BamTools {
namespace Constants {

const uint8_t BAM_SIZEOF_INT = 4;

// header magic number
const char* const BAM_HEADER_MAGIC = "BAM\1";
const uint8_t BAM_HEADER_MAGIC_LENGTH = 4;

// BAM alignment core size
const uint8_t BAM_CORE_SIZE        = 32;
const uint8_t BAM_CORE_BUFFER_SIZE = 8;

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
const uint8_t BAM_CIGAR_MATCH    = 0;
const uint8_t BAM_CIGAR_INS      = 1;
const uint8_t BAM_CIGAR_DEL      = 2;
const uint8_t BAM_CIGAR_REFSKIP  = 3;
const uint8_t BAM_CIGAR_SOFTCLIP = 4;
const uint8_t BAM_CIGAR_HARDCLIP = 5;
const uint8_t BAM_CIGAR_PAD      = 6;
const uint8_t BAM_CIGAR_SEQMATCH = 7;
const uint8_t BAM_CIGAR_MISMATCH = 8;

const char BAM_CIGAR_MATCH_CHAR    = 'M';
const char BAM_CIGAR_INS_CHAR      = 'I';
const char BAM_CIGAR_DEL_CHAR      = 'D';
const char BAM_CIGAR_REFSKIP_CHAR  = 'N';
const char BAM_CIGAR_SOFTCLIP_CHAR = 'S';
const char BAM_CIGAR_HARDCLIP_CHAR = 'H';
const char BAM_CIGAR_PAD_CHAR      = 'P';
const char BAM_CIGAR_SEQMATCH_CHAR = '=';
const char BAM_CIGAR_MISMATCH_CHAR = 'X';

const int BAM_CIGAR_SHIFT = 4;
const int BAM_CIGAR_MASK  = ((1 << BAM_CIGAR_SHIFT) - 1);

// BAM tag types & sizes
const char BAM_TAG_TYPE_ASCII  = 'A';
const char BAM_TAG_TYPE_INT8   = 'c';
const char BAM_TAG_TYPE_UINT8  = 'C';
const char BAM_TAG_TYPE_INT16  = 's';
const char BAM_TAG_TYPE_UINT16 = 'S';
const char BAM_TAG_TYPE_INT32  = 'i';
const char BAM_TAG_TYPE_UINT32 = 'I';
const char BAM_TAG_TYPE_FLOAT  = 'f';
const char BAM_TAG_TYPE_STRING = 'Z';
const char BAM_TAG_TYPE_HEX    = 'H';
const char BAM_TAG_TYPE_ARRAY  = 'B';

const uint8_t BAM_TAG_TAGSIZE        = 2;
const uint8_t BAM_TAG_TYPESIZE       = 1;
const uint8_t BAM_TAG_ARRAYBASE_SIZE = 8;

// DNA bases
const char* const BAM_DNA_LOOKUP = "=ACMGRSVTWYHKDBN";
const uint8_t BAM_BASECODE_EQUAL = 0;
const uint8_t BAM_BASECODE_A     = 1;
const uint8_t BAM_BASECODE_C     = 2;
const uint8_t BAM_BASECODE_M     = 3;
const uint8_t BAM_BASECODE_G     = 4;
const uint8_t BAM_BASECODE_R     = 5;
const uint8_t BAM_BASECODE_S     = 6;
const uint8_t BAM_BASECODE_V     = 7;
const uint8_t BAM_BASECODE_T     = 8;
const uint8_t BAM_BASECODE_W     = 9;
const uint8_t BAM_BASECODE_Y     = 10;
const uint8_t BAM_BASECODE_H     = 11;
const uint8_t BAM_BASECODE_K     = 12;
const uint8_t BAM_BASECODE_D     = 13;
const uint8_t BAM_BASECODE_B     = 14;
const uint8_t BAM_BASECODE_N     = 15;

const char BAM_DNA_EQUAL = '=';
const char BAM_DNA_A     = 'A';
const char BAM_DNA_C     = 'C';
const char BAM_DNA_M     = 'M';
const char BAM_DNA_G     = 'G';
const char BAM_DNA_R     = 'R';
const char BAM_DNA_S     = 'S';
const char BAM_DNA_V     = 'V';
const char BAM_DNA_T     = 'T';
const char BAM_DNA_W     = 'W';
const char BAM_DNA_Y     = 'Y';
const char BAM_DNA_H     = 'H';
const char BAM_DNA_K     = 'K';
const char BAM_DNA_D     = 'D';
const char BAM_DNA_B     = 'B';
const char BAM_DNA_N     = 'N';
const char BAM_DNA_DEL   = '-';
const char BAM_DNA_PAD   = '*';

// zlib & BGZF constants
const char GZIP_ID1   = 31;
const char GZIP_ID2   = 139;
const char CM_DEFLATE = 8;
const char FLG_FEXTRA = 4;
const char OS_UNKNOWN = 255;
const char BGZF_XLEN  = 6;
const char BGZF_ID1   = 66;
const char BGZF_ID2   = 67;
const char BGZF_LEN   = 2;

const int8_t   GZIP_WINDOW_BITS          = -15;
const int8_t   Z_DEFAULT_MEM_LEVEL       = 8;
const uint8_t  BGZF_BLOCK_HEADER_LENGTH  = 18;
const uint8_t  BGZF_BLOCK_FOOTER_LENGTH  = 8;
const uint32_t BGZF_MAX_BLOCK_SIZE       = 65536;
const uint32_t BGZF_DEFAULT_BLOCK_SIZE   = 65536;

} // namespace Constants

//! \cond
// -------------------------
// tag-type helper structs
// -------------------------

// fail on any types not specified below
template<typename T>
struct TagTypeHelper {
    static bool CanConvertFrom(const char) { assert(false); return false; }
    static bool CanConvertTo(const char) { assert(false); return false; }
    static char TypeCode(void) { assert(false); return 0; }
};

template<>
struct TagTypeHelper<uint8_t> {
    static bool CanConvertFrom(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII ||
                 c == Constants::BAM_TAG_TYPE_UINT8 );
    }
    static bool CanConvertTo(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII  ||
                 c == Constants::BAM_TAG_TYPE_UINT8  ||
                 c == Constants::BAM_TAG_TYPE_UINT16 ||
                 c == Constants::BAM_TAG_TYPE_UINT32 );
    }

    static char TypeCode(void) { return Constants::BAM_TAG_TYPE_UINT8; }
};

template<>
struct TagTypeHelper<int8_t> {
    static bool CanConvertFrom(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII ||
                 c == Constants::BAM_TAG_TYPE_INT8 );
    }
    static bool CanConvertTo(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII ||
                 c == Constants::BAM_TAG_TYPE_INT8  ||
                 c == Constants::BAM_TAG_TYPE_INT16 ||
                 c == Constants::BAM_TAG_TYPE_INT32 );
    }
    static char TypeCode(void) { return Constants::BAM_TAG_TYPE_INT8; }
};

template<>
struct TagTypeHelper<uint16_t> {
    static bool CanConvertFrom(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII ||
                 c == Constants::BAM_TAG_TYPE_UINT8 ||
                 c == Constants::BAM_TAG_TYPE_UINT16 );
    }
    static bool CanConvertTo(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_UINT16 ||
                 c == Constants::BAM_TAG_TYPE_UINT32);
    }
    static char TypeCode(void) { return Constants::BAM_TAG_TYPE_UINT16; }
};

template<>
struct TagTypeHelper<int16_t> {
    static bool CanConvertFrom(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII ||
                 c == Constants::BAM_TAG_TYPE_INT8 ||
                 c == Constants::BAM_TAG_TYPE_INT16 );
    }
    static bool CanConvertTo(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_INT16 ||
                 c == Constants::BAM_TAG_TYPE_INT32);
    }
    static char TypeCode(void) { return Constants::BAM_TAG_TYPE_INT16; }
};

template<>
struct TagTypeHelper<uint32_t> {
    static bool CanConvertFrom(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII  ||
                 c == Constants::BAM_TAG_TYPE_UINT8  ||
                 c == Constants::BAM_TAG_TYPE_UINT16 ||
                 c == Constants::BAM_TAG_TYPE_UINT32 );
    }
    static bool CanConvertTo(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_UINT32 );
    }
    static char TypeCode(void) { return Constants::BAM_TAG_TYPE_UINT32; }
};

template<>
struct TagTypeHelper<int32_t> {
    static bool CanConvertFrom(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII  ||
                 c == Constants::BAM_TAG_TYPE_INT8  ||
                 c == Constants::BAM_TAG_TYPE_INT16 ||
                 c == Constants::BAM_TAG_TYPE_INT32 );
    }
    static bool CanConvertTo(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_INT32 );
    }
    static char TypeCode(void) { return Constants::BAM_TAG_TYPE_INT32; }
};

template<>
struct TagTypeHelper<float> {
    static bool CanConvertFrom(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_ASCII  ||
                 c == Constants::BAM_TAG_TYPE_UINT8  ||
                 c == Constants::BAM_TAG_TYPE_INT8   ||
                 c == Constants::BAM_TAG_TYPE_UINT16 ||
                 c == Constants::BAM_TAG_TYPE_INT16  ||
                 c == Constants::BAM_TAG_TYPE_UINT32 ||
                 c == Constants::BAM_TAG_TYPE_INT32  ||
                 c == Constants::BAM_TAG_TYPE_FLOAT);
    }
    static bool CanConvertTo(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_FLOAT );
    }
    static char TypeCode(void) { return Constants::BAM_TAG_TYPE_FLOAT; }
};

template<>
struct TagTypeHelper<std::string> {
    static bool CanConvertFrom(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_HEX ||
                 c == Constants::BAM_TAG_TYPE_STRING );
    }
    static bool CanConvertTo(const char c) {
        return ( c == Constants::BAM_TAG_TYPE_HEX ||
                 c == Constants::BAM_TAG_TYPE_STRING );
    }
    static char TypeCode(void) { return Constants::BAM_TAG_TYPE_STRING; }
};

//! \endcond

} // namespace BamTools

#endif // BAM_CONSTANTS_H
