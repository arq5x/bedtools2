/*  hts_defs.h -- Miscellaneous definitions.

    Copyright (C) 2013-2015,2017 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

#ifndef HTSLIB_HTS_DEFS_H
#define HTSLIB_HTS_DEFS_H

#ifdef __clang__
#ifdef __has_attribute
#define HTS_COMPILER_HAS(attribute) __has_attribute(attribute)
#endif

#elif defined __GNUC__
#define HTS_GCC_AT_LEAST(major, minor) \
    (__GNUC__ > (major) || (__GNUC__ == (major) && __GNUC_MINOR__ >= (minor)))
#endif

#ifndef HTS_COMPILER_HAS
#define HTS_COMPILER_HAS(attribute) 0
#endif
#ifndef HTS_GCC_AT_LEAST
#define HTS_GCC_AT_LEAST(major, minor) 0
#endif

#if HTS_COMPILER_HAS(__noreturn__) || HTS_GCC_AT_LEAST(3,0)
#define HTS_NORETURN __attribute__ ((__noreturn__))
#else
#define HTS_NORETURN
#endif

// GCC introduced warn_unused_result in 3.4 but added -Wno-unused-result later
#if HTS_COMPILER_HAS(__warn_unused_result__) || HTS_GCC_AT_LEAST(4,5)
#define HTS_RESULT_USED __attribute__ ((__warn_unused_result__))
#else
#define HTS_RESULT_USED
#endif

#if HTS_COMPILER_HAS(__unused__) || HTS_GCC_AT_LEAST(3,0)
#define HTS_UNUSED __attribute__ ((__unused__))
#else
#define HTS_UNUSED
#endif

#if HTS_COMPILER_HAS(__deprecated__) || HTS_GCC_AT_LEAST(4,5)
#define HTS_DEPRECATED(message) __attribute__ ((__deprecated__ (message)))
#elif HTS_GCC_AT_LEAST(3,1)
#define HTS_DEPRECATED(message) __attribute__ ((__deprecated__))
#else
#define HTS_DEPRECATED(message)
#endif

#if HTS_COMPILER_HAS(__deprecated__) || HTS_GCC_AT_LEAST(6,4)
#define HTS_DEPRECATED_ENUM(message) __attribute__ ((__deprecated__ (message)))
#else
#define HTS_DEPRECATED_ENUM(message)
#endif

// On mingw the "printf" format type doesn't work.  It needs "gnu_printf"
// in order to check %lld and %z, otherwise it defaults to checking against
// the Microsoft library printf format options despite linking against the
// GNU posix implementation of printf.  The __MINGW_PRINTF_FORMAT macro
// expands to printf or gnu_printf as required, but obviously may not
// exist
#ifdef __MINGW_PRINTF_FORMAT
#define HTS_PRINTF_FMT __MINGW_PRINTF_FORMAT
#else
#define HTS_PRINTF_FMT printf
#endif

#if HTS_COMPILER_HAS(__format__) || HTS_GCC_AT_LEAST(3,0)
#define HTS_FORMAT(type, idx, first) __attribute__((__format__ (type, idx, first)))
#else
#define HTS_FORMAT(type, idx, first)
#endif

#endif
