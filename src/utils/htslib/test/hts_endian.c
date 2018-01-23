/*  test/hts_endian.c -- hts_endian.h unit tests

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "htslib/hts_endian.h"

typedef struct {
    uint8_t u8[2];
    uint8_t u8_unaligned[3];
    int16_t  i16;
    uint16_t u16;
} Test16;

typedef struct {
    uint8_t u8[4];
    uint8_t u8_unaligned[5];
    int32_t  i32;
    uint32_t u32;
} Test32;

typedef struct {
    uint8_t u8[8];
    uint8_t u8_unaligned[9];
    int64_t  i64;
    uint64_t u64;
} Test64;

typedef struct {
    uint8_t u8[4];
    uint8_t u8_unaligned[5];
    float f;
} Test_float;

typedef struct {
    uint8_t u8[8];
    uint8_t u8_unaligned[9];
    double d;
} Test_double;

#define T16(b0, b1, sgn, unsgn) { { b0, b1 }, { 0x00, b0, b1 }, sgn, unsgn }

Test16 tests_16_bit[] = {
    T16(0x00, 0x00,      0,     0),
    T16(0x01, 0x00,      1,     1),
    T16(0x00, 0x01,    256,   256),
    T16(0xff, 0x7f,  32767, 32767),
    T16(0x00, 0x80, -32768, 32768),
    T16(0xff, 0xff,     -1, 65535),
};

#define T32(b0, b1, b2, b3, sgn, unsgn) { \
     { b0, b1, b2, b3 },                  \
     { 0x00, b0, b1, b2, b3 },            \
     sgn, unsgn                           \
}

Test32 tests_32_bit[] = {
    T32(0x00, 0x00, 0x00, 0x00,           0,              0),
    T32(0x01, 0x00, 0x00, 0x00,           1,              1),
    T32(0x00, 0x01, 0x00, 0x00,         256,            256),
    T32(0x00, 0x00, 0x01, 0x00,       65536,          65536),
    T32(0xff, 0xff, 0xff, 0x7f,  2147483647,     2147483647),
    // Odd coding of signed result below avoids a compiler warning
    // as 2147483648 can't fit in a signed 32-bit number
    T32(0x00, 0x00, 0x00, 0x80, -2147483647 - 1, 2147483648U),
    T32(0xff, 0xff, 0xff, 0xff,          -1,     4294967295U),
};

#define T64(b0, b1, b2, b3, b4, b5, b6, b7, sgn, unsgn) { \
     { b0, b1, b2, b3, b4, b5, b6, b7 },                  \
     { 0x00, b0, b1, b2, b3, b4, b5, b6, b7 },            \
     sgn, unsgn                                           \
}


Test64 tests_64_bit[] = {
    T64(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0, 0),
    T64(0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 1, 1),
    T64(0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 256, 256),
    T64(0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 65536, 65536),
    T64(0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 4294967296, 4294967296),
    T64(0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f,
        9223372036854775807LL, 9223372036854775807ULL),
    // Odd coding of signed result below avoids a compiler warning
    // as 9223372036854775808 can't fit in a signed 64-bit number
    T64(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80,
        -9223372036854775807LL - 1LL, 9223372036854775808ULL),
    T64(0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        -1, 18446744073709551615ULL),
};

#define TF(b0, b1, b2, b3, f) { { b0, b1, b2, b3 }, { 0x00, b0, b1, b2, b3}, f }

Test_float tests_float[] = {
    TF(0x00, 0x00, 0x00, 0x00,  0.0f),
    TF(0x00, 0x00, 0x80, 0x3f,  1.0f),
    TF(0x00, 0x00, 0x80, 0xbf, -1.0f),
    TF(0x00, 0x00, 0x20, 0x41, 10.0f),
    TF(0xd0, 0x0f, 0x49, 0x40,  3.14159f),
    TF(0xa8, 0x0a, 0xff, 0x66,  6.022e23f),
    TF(0xcd, 0x84, 0x03, 0x13,  1.66e-27f),
};

#define TD(b0, b1, b2, b3, b4, b5, b6, b7, d) { \
    { b0, b1, b2, b3, b4, b5, b6, b7 },         \
    { 0x00, b0, b1, b2, b3, b4, b5, b6, b7 },   \
    d                                           \
}

Test_double tests_double[] = {
    TD(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  0.0),
    TD(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,  1.0),
    TD(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0xbf, -1.0),
    TD(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x40, 10.0),
    TD(0x18, 0x2d, 0x44, 0x54, 0xfb, 0x21, 0x09, 0x40,  3.141592653589793),
    TD(0x2b, 0x08, 0x0c, 0xd3, 0x85, 0xe1, 0xdf, 0x44,  6.022140858e23),
    TD(0x55, 0xfa, 0x81, 0x74, 0xf7, 0x71, 0x60, 0x3a,  1.66053904e-27),
};

#define NELE(x) (sizeof(x)/sizeof(x[0]))

static char * to_hex(uint8_t *buf, int len) {
    static char str[64];
    int i, o;

    for (i = 0, o = 0; i < len; i++, o += 3) {
        snprintf(str + o, sizeof(str) - o, "%02x ", buf[i]);
    }
    return str;
}

static int t16_bit(int verbose) {
    uint8_t buf[9];
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_16_bit); i++) {
        uint16_t u16;
        int16_t  i16;

        if (verbose) {
            fprintf(stderr, "%s %6"PRId16" %6"PRId16"\n",
                    to_hex(tests_16_bit[i].u8, 2),
                    tests_16_bit[i].i16, tests_16_bit[i].u16);
        }

        u16 = le_to_u16(tests_16_bit[i].u8);
        if (u16 != tests_16_bit[i].u16) {
            fprintf(stderr, "Failed %s => %"PRIu16"; expected %"PRIu16"\n",
                    to_hex(tests_16_bit[i].u8, 2), u16, tests_16_bit[i].u16);
            errors++;
        }

        i16 = le_to_i16(tests_16_bit[i].u8);
        if (i16 != tests_16_bit[i].i16) {
            fprintf(stderr, "Failed %s => %"PRId16"; expected %"PRId16"\n",
                    to_hex(tests_16_bit[i].u8, 2), i16, tests_16_bit[i].i16);
            errors++;
        }

        u16 = le_to_u16(tests_16_bit[i].u8_unaligned + 1);
        if (u16 != tests_16_bit[i].u16) {
            fprintf(stderr,
                    "Failed unaligned %s => %"PRIu16"; expected %"PRIu16"\n",
                    to_hex(tests_16_bit[i].u8_unaligned + 1, 2),
                    u16, tests_16_bit[i].u16);
            errors++;
        }

        i16 = le_to_i16(tests_16_bit[i].u8_unaligned + 1);
        if (i16 != tests_16_bit[i].i16) {
            fprintf(stderr,
                    "Failed unaligned %s => %"PRId16"; expected %"PRId16"\n",
                    to_hex(tests_16_bit[i].u8_unaligned + 1, 2),
                    i16, tests_16_bit[i].i16);
            errors++;
        }

        u16_to_le(tests_16_bit[i].u16, buf);
        if (memcmp(buf, tests_16_bit[i].u8, 2) != 0) {
            fprintf(stderr, "Failed %"PRIu16" => %s; expected %s\n",
                    tests_16_bit[i].u16,
                    to_hex(buf, 2), to_hex(tests_16_bit[i].u8, 2));
            errors++;
        }

        i16_to_le(tests_16_bit[i].i16, buf);
        if (memcmp(buf, tests_16_bit[i].u8, 2) != 0) {
            fprintf(stderr, "Failed %"PRId16" => %s; expected %s\n",
                    tests_16_bit[i].i16,
                    to_hex(buf, 2), to_hex(tests_16_bit[i].u8, 2));
            errors++;
        }

        u16_to_le(tests_16_bit[i].u16, buf + 1);
        if (memcmp(buf + 1, tests_16_bit[i].u8, 2) != 0) {
            fprintf(stderr, "Failed unaligned %"PRIu16" => %s; expected %s\n",
                    tests_16_bit[i].u16,
                    to_hex(buf + 1, 2), to_hex(tests_16_bit[i].u8, 2));
            errors++;
        }

        i16_to_le(tests_16_bit[i].i16, buf + 1);
        if (memcmp(buf + 1, tests_16_bit[i].u8, 2) != 0) {
            fprintf(stderr, "Failed unaligned %"PRId16" => %s; expected %s\n",
                    tests_16_bit[i].i16,
                    to_hex(buf + 1, 2), to_hex(tests_16_bit[i].u8, 2));
            errors++;
        }
    }

    return errors;
}

static int t32_bit(int verbose) {
    uint8_t buf[9];
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_32_bit); i++) {
        uint32_t u32;
        int32_t  i32;

        if (verbose) {
            fprintf(stderr, "%s %11"PRId32" %11"PRIu32"\n",
                    to_hex(tests_32_bit[i].u8, 4),
                    tests_32_bit[i].i32, tests_32_bit[i].u32);
        }

        u32 = le_to_u32(tests_32_bit[i].u8);
        if (u32 != tests_32_bit[i].u32) {
            fprintf(stderr, "Failed %s => %"PRIu32"; expected %"PRIu32"\n",
                    to_hex(tests_32_bit[i].u8, 4), u32, tests_32_bit[i].u32);
            errors++;
        }
        i32 = le_to_i32(tests_32_bit[i].u8);
        if (i32 != tests_32_bit[i].i32) {
            fprintf(stderr, "Failed %s => %"PRId32"; expected %"PRId32"\n",
                    to_hex(tests_32_bit[i].u8, 4), i32, tests_32_bit[i].i32);
            errors++;
        }

        u32 = le_to_u32(tests_32_bit[i].u8_unaligned + 1);
        if (u32 != tests_32_bit[i].u32) {
            fprintf(stderr,
                    "Failed unaligned %s => %"PRIu32"; expected %"PRIu32"\n",
                    to_hex(tests_32_bit[i].u8_unaligned + 1, 4),
                    u32, tests_32_bit[i].u32);
            errors++;
        }
        i32 = le_to_i32(tests_32_bit[i].u8_unaligned + 1);
        if (i32 != tests_32_bit[i].i32) {
            fprintf(stderr,
                    "Failed unaligned %s => %"PRId32"; expected %"PRId32"\n",
                    to_hex(tests_32_bit[i].u8_unaligned + 1, 4),
                    i32, tests_32_bit[i].i32);
            errors++;
        }

        u32_to_le(tests_32_bit[i].u32, buf);
        if (memcmp(buf, tests_32_bit[i].u8, 4) != 0) {
            fprintf(stderr, "Failed %"PRIu32" => %s; expected %s\n",
                    tests_32_bit[i].u32,
                    to_hex(buf, 4), to_hex(tests_32_bit[i].u8, 4));
            errors++;
        }

        i32_to_le(tests_32_bit[i].i32, buf);
        if (memcmp(buf, tests_32_bit[i].u8, 4) != 0) {
            fprintf(stderr, "Failed %"PRId32" => %s; expected %s\n",
                    tests_32_bit[i].i32,
                    to_hex(buf, 4), to_hex(tests_32_bit[i].u8, 4));
            errors++;
        }

        u32_to_le(tests_32_bit[i].u32, buf + 1);
        if (memcmp(buf + 1, tests_32_bit[i].u8, 4) != 0) {
            fprintf(stderr, "Failed unaligned %"PRIu32" => %s; expected %s\n",
                    tests_32_bit[i].u32,
                    to_hex(buf + 1, 4), to_hex(tests_32_bit[i].u8, 4));
            errors++;
        }

        i32_to_le(tests_32_bit[i].i32, buf + 1);
        if (memcmp(buf + 1, tests_32_bit[i].u8, 4) != 0) {
            fprintf(stderr, "Failed unaligned %"PRId32" => %s; expected %s\n",
                    tests_32_bit[i].i32,
                    to_hex(buf + 1, 4), to_hex(tests_32_bit[i].u8, 4));
            errors++;
        }
    }

    return errors;
}

static int t64_bit(int verbose) {
    uint8_t buf[9];
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_64_bit); i++) {
        uint64_t u64;
        int64_t  i64;

        if (verbose) {
            fprintf(stderr, "%s %20"PRId64" %20"PRIu64"\n",
                    to_hex(tests_64_bit[i].u8, 8),
                    tests_64_bit[i].i64, tests_64_bit[i].u64);
        }

        u64 = le_to_u64(tests_64_bit[i].u8);
        if (u64 != tests_64_bit[i].u64) {
            fprintf(stderr, "Failed %s => %"PRIu64"; expected %"PRIu64"\n",
                    to_hex(tests_64_bit[i].u8, 8), u64, tests_64_bit[i].u64);
            errors++;
        }

        i64 = le_to_i64(tests_64_bit[i].u8);
        if (i64 != tests_64_bit[i].i64) {
            fprintf(stderr, "Failed %s => %"PRId64"; expected %"PRId64"\n",
                    to_hex(tests_64_bit[i].u8, 8), i64, tests_64_bit[i].i64);
            errors++;
        }

        u64 = le_to_u64(tests_64_bit[i].u8_unaligned + 1);
        if (u64 != tests_64_bit[i].u64) {
            fprintf(stderr,
                    "Failed unaligned %s => %"PRIu64"; expected %"PRIu64"\n",
                    to_hex(tests_64_bit[i].u8_unaligned + 1, 8),
                    u64, tests_64_bit[i].u64);
            errors++;
        }

        i64 = le_to_i64(tests_64_bit[i].u8_unaligned + 1);
        if (i64 != tests_64_bit[i].i64) {
            fprintf(stderr,
                    "Failed unaligned %s => %"PRId64"; expected %"PRId64"\n",
                    to_hex(tests_64_bit[i].u8_unaligned + 1, 8),
                    i64, tests_64_bit[i].i64);
            errors++;
        }

        u64_to_le(tests_64_bit[i].u64, buf);
        if (memcmp(buf, tests_64_bit[i].u8, 8) != 0) {
            fprintf(stderr, "Failed %"PRIu64" => %s; expected %s\n",
                    tests_64_bit[i].u64,
                    to_hex(buf, 8), to_hex(tests_64_bit[i].u8, 8));
            errors++;
        }

        i64_to_le(tests_64_bit[i].i64, buf);
        if (memcmp(buf, tests_64_bit[i].u8, 8) != 0) {
            fprintf(stderr, "Failed %"PRId64" => %s; expected %s\n",
                    tests_64_bit[i].i64,
                    to_hex(buf, 8), to_hex(tests_64_bit[i].u8, 8));
            errors++;
        }

        u64_to_le(tests_64_bit[i].u64, buf + 1);
        if (memcmp(buf + 1, tests_64_bit[i].u8, 8) != 0) {
            fprintf(stderr, "Failed unaligned %"PRIu64" => %s; expected %s\n",
                    tests_64_bit[i].u64,
                    to_hex(buf + 1, 8), to_hex(tests_64_bit[i].u8, 8));
            errors++;
        }

        i64_to_le(tests_64_bit[i].i64, buf + 1);
        if (memcmp(buf + 1, tests_64_bit[i].u8, 8) != 0) {
            fprintf(stderr, "Failed unaligned %"PRId64" => %s; expected %s\n",
                    tests_64_bit[i].i64,
                    to_hex(buf + 1, 8), to_hex(tests_64_bit[i].u8, 8));
            errors++;
        }
    }

    return errors;
}

int t_float(int verbose) {
    uint8_t buf[9];
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_float); i++) {
        float f;

        if (verbose) {
            fprintf(stderr, "%s %g\n",
                    to_hex(tests_float[i].u8, 4), tests_float[i].f);
        }

        f = le_to_float(tests_float[i].u8);
        if (f != tests_float[i].f) {
            fprintf(stderr, "Failed %s => %g; expected %g\n",
                    to_hex(tests_float[i].u8, 4), f, tests_float[i].f);
            errors++;
        }

        f = le_to_float(tests_float[i].u8_unaligned + 1);
        if (f != tests_float[i].f) {
            fprintf(stderr, "Failed unaligned %s => %g; expected %g\n",
                    to_hex(tests_float[i].u8_unaligned + 1, 4),
                    f, tests_float[i].f);
            errors++;
        }

        float_to_le(tests_float[i].f, buf);
        if (memcmp(tests_float[i].u8, buf, 4) != 0) {
            fprintf(stderr, "Failed %g => %s; expected %s\n",
                    tests_float[i].f,
                    to_hex(buf, 4), to_hex(tests_float[i].u8, 4));
        }

        float_to_le(tests_float[i].f, buf + 1);
        if (memcmp(tests_float[i].u8, buf + 1, 4) != 0) {
            fprintf(stderr, "Failed unaligned %g => %s; expected %s\n",
                    tests_float[i].f,
                    to_hex(buf + 1, 4), to_hex(tests_float[i].u8, 4));
        }
    }
    return errors;
}

int t_double(int verbose) {
    uint8_t buf[9];
    size_t i;
    int errors = 0;

    for (i = 0; i < NELE(tests_double); i++) {
        double f;

        if (verbose) {
            fprintf(stderr, "%s %.15g\n",
                    to_hex(tests_double[i].u8, 8), tests_double[i].d);
        }

        f = le_to_double(tests_double[i].u8);
        if (f != tests_double[i].d) {
            fprintf(stderr, "Failed %s => %.15g; expected %.15g\n",
                    to_hex(tests_double[i].u8, 8), f, tests_double[i].d);
            errors++;
        }

        f = le_to_double(tests_double[i].u8_unaligned + 1);
        if (f != tests_double[i].d) {
            fprintf(stderr, "Failed unaligned %s => %.15g; expected %.15g\n",
                    to_hex(tests_double[i].u8_unaligned + 1, 8),
                    f, tests_double[i].d);
            errors++;
        }

        double_to_le(tests_double[i].d, buf);
        if (memcmp(tests_double[i].u8, buf, 8) != 0) {
            fprintf(stderr, "Failed %.15g => %s; expected %s\n",
                    tests_double[i].d,
                    to_hex(buf, 8), to_hex(tests_double[i].u8, 8));
        }

        double_to_le(tests_double[i].d, buf + 1);
        if (memcmp(tests_double[i].u8, buf + 1, 8) != 0) {
            fprintf(stderr, "Failed unaligned %.15g => %s; expected %s\n",
                    tests_double[i].d,
                    to_hex(buf + 1, 8), to_hex(tests_double[i].u8, 8));
        }
    }
    return errors;
}

int main(int argc, char **argv) {
    int verbose = 0;
    int errors = 0;

    if (argc > 1 && strcmp(argv[1], "-v") == 0) verbose = 1;

    errors += t16_bit(verbose);
    errors += t32_bit(verbose);
    errors += t64_bit(verbose);
    errors += t_float(verbose);
    errors += t_double(verbose);
    if (errors) {
        fprintf(stderr, "%d errors\n", errors);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
