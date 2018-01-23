/// @file hts_endian.h
/// Byte swapping and unaligned access functions.
/*
   Copyright (C) 2017 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTS_ENDIAN_H
#define HTS_ENDIAN_H

#include <stdint.h>

/*
 * Compile-time endianness tests.
 *
 * Note that these tests may fail.  They should only be used to enable
 * faster versions of endian-neutral implementations.  The endian-neutral
 * version should always be available as a fall-back.
 *
 * See https://sourceforge.net/p/predef/wiki/Endianness/
 */

/* Save typing as both endian and unaligned tests want to know about x86 */
#if (defined(__i386__) || defined(__i386) || defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64) || defined(__i686__) || defined(__i686)) && !defined(HTS_x86)
#    define HTS_x86  /* x86 and x86_64 platform */
#endif

/** @def HTS_LITTLE_ENDIAN
 *  @brief Defined if platform is known to be little-endian
 */

#ifndef HTS_LITTLE_ENDIAN
#    if (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) \
      || defined(__LITTLE_ENDIAN__) \
      || defined(HTS_x86) \
      || defined(__ARMEL__) || defined(__THUMBEL__) || defined(__AARCH64EL__) \
      || defined(_MIPSEL) || defined(__MIPSEL) || defined(__MIPSEL__)
#        define HTS_LITTLE_ENDIAN
#    endif
#endif

/** @def HTS_BIG_ENDIAN
 *  @brief Defined if platform is known to be big-endian
 */

#ifndef HTS_BIG_ENDIAN
#    if (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__) \
      || defined(__BIG_ENDIAN__) \
      || defined(__ARMEB__) || defined(__THUMBEB__) || defined(__AAARCHEB__) \
      || defined(_MIPSEB) || defined(__MIPSEB) || defined(__MIPSEB__)
#        define HTS_BIG_ENDIAN
#    endif
#endif

/** @def HTS_ENDIAN_NEUTRAL
 *  @brief Define this to disable any endian-specific optimizations
 */

#if defined(HTS_ENDIAN_NEUTRAL) || (defined(HTS_LITTLE_ENDIAN) && defined(HTS_BIG_ENDIAN))
/* Disable all endian-specific code. */
#    undef HTS_LITTLE_ENDIAN
#    undef HTS_BIG_ENDIAN
#endif

/** @def HTS_ALLOW_UNALIGNED
 *  @brief Control use of unaligned memory access.
 *
 * Defining HTS_ALLOW_UNALIGNED=1 converts shift-and-or to simple casts on
 * little-endian platforms that can tolerate unaligned access (notably Intel
 * x86).
 *
 * Defining HTS_ALLOW_UNALIGNED=0 forces shift-and-or.
 */

// Consider using AX_CHECK_ALIGNED_ACCESS_REQUIRED in autoconf.
#ifndef HTS_ALLOW_UNALIGNED
#    if defined(HTS_x86)
#        define HTS_ALLOW_UNALIGNED 1
#    else
#        define HTS_ALLOW_UNALIGNED 0
#    endif
#endif

#if HTS_ALLOW_UNALIGNED != 0
#    if defined (__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3))
// This prevents problems with gcc's vectoriser generating the wrong
// instructions for unaligned data.
typedef uint16_t uint16_u __attribute__ ((__aligned__ (1)));
typedef uint32_t uint32_u __attribute__ ((__aligned__ (1)));
typedef uint64_t uint64_u __attribute__ ((__aligned__ (1)));
#else
typedef uint16_t uint16_u;
typedef uint32_t uint32_u;
typedef uint64_t uint64_u;
#    endif
#endif

/// Get a uint16_t value from an unsigned byte array
/** @param buf Pointer to source byte, may be unaligned
 *  @return A 16 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
static inline uint16_t le_to_u16(const uint8_t *buf) {
#if defined(HTS_LITTLE_ENDIAN) && HTS_ALLOW_UNALIGNED != 0
    return *((uint16_u *) buf);
#else
    return (uint16_t) buf[0] | ((uint16_t) buf[1] << 8);
#endif
}

/// Get a uint32_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
static inline uint32_t le_to_u32(const uint8_t *buf) {
#if defined(HTS_LITTLE_ENDIAN) && HTS_ALLOW_UNALIGNED != 0
    return *((uint32_u *) buf);
#else
    return ((uint32_t) buf[0] |
            ((uint32_t) buf[1] << 8) |
            ((uint32_t) buf[2] << 16) |
            ((uint32_t) buf[3] << 24));
#endif
}

/// Get a uint64_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit unsigned integer
 *  The input is read in little-endian byte order.
 */
static inline uint64_t le_to_u64(const uint8_t *buf) {
#if defined(HTS_LITTLE_ENDIAN) && HTS_ALLOW_UNALIGNED != 0
    return *((uint64_u *) buf);
#else
    return ((uint64_t) buf[0] |
            ((uint64_t) buf[1] << 8) |
            ((uint64_t) buf[2] << 16) |
            ((uint64_t) buf[3] << 24) |
            ((uint64_t) buf[4] << 32) |
            ((uint64_t) buf[5] << 40) |
            ((uint64_t) buf[6] << 48) |
            ((uint64_t) buf[7] << 56));
#endif
}

/// Store a uint16_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
static inline void u16_to_le(uint16_t val, uint8_t *buf) {
#if defined(HTS_LITTLE_ENDIAN) && HTS_ALLOW_UNALIGNED != 0
    *((uint16_u *) buf) = val;
#else
    buf[0] = val & 0xff;
    buf[1] = (val >> 8) & 0xff;
#endif
}

/// Store a uint32_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
static inline void u32_to_le(uint32_t val, uint8_t *buf) {
#if defined(HTS_LITTLE_ENDIAN) && HTS_ALLOW_UNALIGNED != 0
    *((uint32_u *) buf) = val;
#else
    buf[0] = val & 0xff;
    buf[1] = (val >> 8) & 0xff;
    buf[2] = (val >> 16) & 0xff;
    buf[3] = (val >> 24) & 0xff;
#endif
}

/// Store a uint64_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
static inline void u64_to_le(uint64_t val, uint8_t *buf) {
#if defined(HTS_LITTLE_ENDIAN) && HTS_ALLOW_UNALIGNED != 0
    *((uint64_u *) buf) = val;
#else
    buf[0] = val & 0xff;
    buf[1] = (val >> 8) & 0xff;
    buf[2] = (val >> 16) & 0xff;
    buf[3] = (val >> 24) & 0xff;
    buf[4] = (val >> 32) & 0xff;
    buf[5] = (val >> 40) & 0xff;
    buf[6] = (val >> 48) & 0xff;
    buf[7] = (val >> 56) & 0xff;
#endif
}

/* Signed values.  Grab the data as unsigned, then convert to signed without
 * triggering undefined behaviour.  On any sensible platform, the conversion
 * should optimise away to nothing.
 */

/// Get an int8_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 8 bit signed integer
 *  The input data is interpreted as 2's complement representation.
 */
static inline int8_t le_to_i8(const uint8_t *buf) {
    return *buf < 0x80 ? *buf : -((int8_t) (0xff - *buf)) - 1;
}

/// Get an int16_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 16 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
static inline int16_t le_to_i16(const uint8_t *buf) {
    uint16_t v = le_to_u16(buf);
    return v < 0x8000 ? v : -((int16_t) (0xffff - v)) - 1;
}

/// Get an int32_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
static inline int32_t le_to_i32(const uint8_t *buf) {
    uint32_t v = le_to_u32(buf);
    return v < 0x80000000U ? v : -((int32_t) (0xffffffffU - v)) - 1;
}

/// Get an int64_t value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit signed integer
 *  The input data is interpreted as 2's complement representation in
 *  little-endian byte order.
 */
static inline int64_t le_to_i64(const uint8_t *buf) {
    uint64_t v = le_to_u64(buf);
    return (v < 0x8000000000000000ULL
            ? v : -((int64_t) (0xffffffffffffffffULL - v)) - 1);
}

// Converting the other way is easier as signed -> unsigned is well defined.

/// Store a uint16_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
static inline void i16_to_le(int16_t val, uint8_t *buf) {
    u16_to_le(val, buf);
}

/// Store a uint32_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
static inline void i32_to_le(int32_t val, uint8_t *buf) {
    u32_to_le(val, buf);
}

/// Store a uint64_t value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
static inline void i64_to_le(int64_t val, uint8_t *buf) {
    u64_to_le(val, buf);
}

/* Floating point.  Assumptions:
 *  Platform uses IEEE 754 format
 *  sizeof(float) == sizeof(uint32_t)
 *  sizeof(double) == sizeof(uint64_t)
 *  Endian-ness is the same for both floating point and integer
 *  Type-punning via a union is allowed
 */

/// Get a float value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 32 bit floating point value
 *  The input is interpreted as an IEEE 754 format float in little-endian
 *  byte order.
 */
static inline float le_to_float(const uint8_t *buf) {
    union {
        uint32_t u;
        float   f;
    } convert;

    convert.u = le_to_u32(buf);
    return convert.f;
}

/// Get a double value from an unsigned byte array
/** @param buf Pointer to source byte array, may be unaligned
 *  @return A 64 bit floating point value
 *  The input is interpreted as an IEEE 754 format double in little-endian
 *  byte order.
 */
static inline double le_to_double(const uint8_t *buf) {
    union {
        uint64_t u;
        double   f;
    } convert;

    convert.u = le_to_u64(buf);
    return convert.f;
}

/// Store a float value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
static inline void float_to_le(float val, uint8_t *buf) {
    union {
        uint32_t u;
        float f;
    } convert;

    convert.f = val;
    u32_to_le(convert.u, buf);
}

/// Store a double value in little-endian byte order
/** @param val The value to store
 *  @param buf Where to store it (may be unaligned)
 */
static inline void double_to_le(double val, uint8_t *buf) {
    union {
        uint64_t u;
        double f;
    } convert;

    convert.f = val;
    u64_to_le(convert.u, buf);
}

#endif /* HTS_ENDIAN_H */
