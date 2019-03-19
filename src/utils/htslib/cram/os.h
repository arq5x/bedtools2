/*
Copyright (c) 1993, 1995-2002 MEDICAL RESEARCH COUNCIL
All rights reserved

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1 Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2 Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3 Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF
MOLECULAR BIOLOGY nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
Copyright (c) 2004, 2006, 2009-2011, 2013, 2017 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * File: os.h
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *         Hills Road
 *         Cambridge CB2 2QH
 *         United Kingdom
 *
 * Description: operating system specific type definitions
 *
 */

#ifndef _OS_H_
#define _OS_H_

#include <limits.h>
#include <inttypes.h>
#include "htslib/hts_endian.h"

#ifdef __cplusplus
extern "C" {
#endif


/*-----------------------------------------------------------------------------
 * Byte swapping macros
 */

/*
 * Our new swap runs at the same speed on Ultrix, but substantially faster
 * (300% for swap_int4, ~50% for swap_int2) on an Alpha (due to the lack of
 * decent 'char' support).
 *
 * They also have the ability to swap in situ (src == dst). Newer code now
 * relies on this so don't change back!
 */
#define iswap_int8(x)                           \
    (((x & 0x00000000000000ffLL) << 56) +       \
     ((x & 0x000000000000ff00LL) << 40) +       \
     ((x & 0x0000000000ff0000LL) << 24) +       \
     ((x & 0x00000000ff000000LL) <<  8) +       \
     ((x & 0x000000ff00000000LL) >>  8) +       \
     ((x & 0x0000ff0000000000LL) >> 24) +       \
     ((x & 0x00ff000000000000LL) >> 40) +       \
     ((x & 0xff00000000000000LL) >> 56))

#define iswap_int4(x)                           \
    (((x & 0x000000ff) << 24) +                 \
     ((x & 0x0000ff00) <<  8) +                 \
     ((x & 0x00ff0000) >>  8) +                 \
     ((x & 0xff000000) >> 24))

#define iswap_int2(x)                           \
    (((x & 0x00ff) << 8) +                      \
     ((x & 0xff00) >> 8))

/*
 * Linux systems may use byteswap.h to get assembly versions of byte-swap
 * on intel systems. This can be as trivial as the bswap opcode, which works
 * out at over 2-times faster than iswap_int4 above.
 */
#if 0
#if defined(__linux__)
#    include <byteswap.h>
#    undef iswap_int8
#    undef iswap_int4
#    undef iswap_int2
#    define iswap_int8 bswap_64
#    define iswap_int4 bswap_32
#    define iswap_int2 bswap_16
#endif
#endif


/*
 * Macros to specify that data read in is of a particular endianness.
 * The macros here swap to the appropriate order for the particular machine
 * running the macro and return the new answer. These may also be used when
 * writing to a file to specify that we wish to write in (eg) big endian
 * format.
 *
 * This leads to efficient code as most of the time these macros are
 * trivial.
 */
#if defined(HTS_BIG_ENDIAN)
#define le_int4(x) iswap_int4((x))
#define le_int2(x) iswap_int2((x))
#elif defined(HTS_LITTLE_ENDIAN)
#define le_int4(x) (x)
#define le_int2(x) (x)
#else
static inline uint32_t le_int4(uint32_t x) {
    return le_to_u32((uint8_t *) &x);
}
static inline uint16_t le_int2(uint16_t x) {
    return le_to_u16((uint8_t *) &x);
}
#endif

/*-----------------------------------------------------------------------------
 * <inttypes.h> definitions, incase they're not present
 */

#ifndef PRId64
#define __PRI64__ "l"
#define PRId64 __PRI64__ "d"
#define PRId32 "d"
#define PRId16 "d"
#define PRId8  "d"
#define PRIu64 __PRI64__ "u"
#define PRIu32 "u"
#define PRIu16 "u"
#define PRIu8  "u"
#endif

/*-----------------------------------------------------------------------------
 * Operating system specifics.
 * These ought to be done by autoconf, but are legacy code.
 */
/*
 * SunOS 4.x
 * Even though we use the ANSI gcc, we make use the the standard SunOS 4.x
 * libraries and include files, which are non-ansi
 */
#if defined(__sun__) && !defined(__svr4__)
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

/*
 * Microsoft Visual C++
 * Windows
 */
#if defined(_MSC_VER)
#define popen _popen
#define pclose _pclose
#define ftruncate(fd,len) _chsize(fd,len)
#endif


/*
 * Microsoft Windows running MinGW
 */
#if defined(__MINGW32__)
#include <io.h>
#define mkdir(filename,mode) mkdir((filename))
#define sysconf(x) 512
#ifndef ftruncate
#  define ftruncate(fd,len) _chsize(fd,len)
#endif
#endif

#ifdef __cplusplus
}
#endif

#endif /*_OS_H_*/
