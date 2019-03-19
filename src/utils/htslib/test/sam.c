/*  test/sam.c -- SAM/BAM/CRAM API test cases.

    Copyright (C) 2014-2017 Genome Research Ltd.

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

#include <config.h>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

// Suppress message for faidx_fetch_nseq(), which we're intentionally testing
#include "htslib/hts_defs.h"
#undef HTS_DEPRECATED
#define HTS_DEPRECATED(message)

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

int status;

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) fail(const char *fmt, ...)
{
    va_list args;

    fprintf(stderr, "Failed: ");
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");

    status = EXIT_FAILURE;
}

uint8_t *check_bam_aux_get(const bam1_t *aln, const char *tag, char type)
{
    uint8_t *p = bam_aux_get(aln, tag);
    if (p) {
        if (*p == type) return p;
        else fail("%s field of type '%c', expected '%c'\n", tag, *p, type);
    }
    else fail("can't find %s field\n", tag);

    return NULL;
}

static void check_int_B_array(bam1_t *aln, char *tag,
                             uint32_t nvals, int64_t *vals) {
    uint8_t *p;
    if ((p = check_bam_aux_get(aln, tag, 'B')) != NULL) {
        uint32_t i;

        if (bam_auxB_len(p) != nvals)
            fail("Wrong length reported for %s field, got %u, expected %u\n",
                 tag, bam_auxB_len(p), nvals);

        for (i = 0; i < nvals; i++) {
            if (bam_auxB2i(p, i) != vals[i]) {
                fail("Wrong value from bam_auxB2i for %s field index %u, "
                     "got %"PRId64" expected %"PRId64"\n",
                     tag, i, bam_auxB2i(p, i), vals[i]);
            }
            if (bam_auxB2f(p, i) != (double) vals[i]) {
                fail("Wrong value from bam_auxB2f for %s field index %u, "
                     "got %f expected %f\n",
                     tag, i, bam_auxB2f(p, i), (double) vals[i]);
            }
        }
    }
}

#define PI 3.141592653589793
#define E  2.718281828459045
#define HELLO "Hello, world!"
#define NEW_HELLO "Yo, dude"
#define BEEF "DEADBEEF"

#define str(x) #x
#define xstr(x) str(x)

#define NELE(x) (sizeof(x)/sizeof(x[0]))

static int test_update_int(bam1_t *aln,
                           const char target_id[2], int64_t target_val,
                           char expected_type,
                           const char next_id[2], int64_t next_val,
                           char next_type) {
    uint8_t *p;

    // Try updating target
    if (bam_aux_update_int(aln, target_id, target_val) < 0) {
        fail("update %.2s tag", target_id);
        return -1;
    }

    // Check it's there and has the right type and value
    p = bam_aux_get(aln, target_id);
    if (!p) {
        fail("find  %.2s tag", target_id);
        return -1;
    }
    if (*p != expected_type || bam_aux2i(p) != target_val) {
        fail("%.2s field is %c:%"PRId64"; expected %c:%"PRId64,
             target_id, *p, bam_aux2i(p), expected_type, target_val);
        return -1;
    }

    // If given, check that the next tag hasn't been clobbered by the
    // update above.
    if (!*next_id) return 0;
    p = bam_aux_get(aln, next_id);
    if (!p) {
        fail("find  %.2s tag after updating %.2s", next_id, target_id);
        return -1;
    }
    if (*p != next_type || bam_aux2i(p) != next_val) {
        fail("after updating %.2s to %"PRId64":"
             " %.2s field is %c:%"PRId64"; expected %c:%"PRId64,
             target_id, target_val,
             next_id, *p, bam_aux2i(p), next_type, next_val);
        return -1;
    }
    return 0;
}

#define CHECK_ARRAY_VALS(T, GET_VAL, FMT1, FMT2) do {                   \
        T * vals = (T *) data;                                          \
    uint32_t i;                                                         \
    for (i = 0; i < nitems; i++) {                                      \
        if (GET_VAL(p, i) != vals[i]) {                                 \
            fail("Wrong value from %s for %.2s field index %u, "        \
                 "got %" FMT1 " expected %" FMT2,                       \
                 xstr(GET_VAL), target_id, i, GET_VAL(p, i), vals[i]);  \
            return -1;                                                  \
        }                                                               \
    }                                                                   \
} while (0)

static int test_update_array(bam1_t *aln, const char target_id[2],
                             uint8_t type, uint32_t nitems, void *data,
                             const char next_id[2], int64_t next_val,
                             char next_type)
{
    uint8_t *p;

    // Try updating target
    if (bam_aux_update_array(aln, target_id, type, nitems, data) < 0) {
        fail("update %2.s tag", target_id);
        return -1;
    }

    // Check values
    p = bam_aux_get(aln, target_id);
    if (!p) {
        fail("find  %.2s tag", target_id);
        return -1;
    }
    switch (type) {
        case 'c':
            CHECK_ARRAY_VALS(int8_t, bam_auxB2i, PRId64, PRId8); break;
        case 'C':
            CHECK_ARRAY_VALS(uint8_t, bam_auxB2i, PRId64, PRIu8); break;
        case 's':
            CHECK_ARRAY_VALS(int16_t, bam_auxB2i, PRId64, PRId16); break;
        case 'S':
            CHECK_ARRAY_VALS(uint16_t, bam_auxB2i, PRId64, PRIu16); break;
        case 'i':
            CHECK_ARRAY_VALS(int32_t, bam_auxB2i, PRId64, PRId32); break;
        case 'I':
            CHECK_ARRAY_VALS(uint32_t, bam_auxB2i, PRId64, PRIu32); break;
        case 'f':
            CHECK_ARRAY_VALS(float, bam_auxB2f, "e", "e"); break;
    }

    // If given, check that the next tag hasn't been clobbered by the
    // update above.
    if (!*next_id) return 0;
    p = bam_aux_get(aln, next_id);
    if (!p) {
        fail("find  %.2s tag after updating %.2s", next_id, target_id);
        return -1;
    }
    if (*p != next_type || bam_aux2i(p) != next_val) {
        fail("after updating %.2s:"
             " %.2s field is %c:%"PRId64"; expected %c:%"PRId64,
             target_id, next_id, *p, bam_aux2i(p), next_type, next_val);
        return -1;
    }

    return 0;
}

static int aux_fields1(void)
{
    static const char sam[] = "data:,"
"@SQ\tSN:one\tLN:1000\n"
"@SQ\tSN:two\tLN:500\n"
"r1\t0\tone\t500\t20\t8M\t*\t0\t0\tATGCATGC\tqqqqqqqq\tXA:A:k\tXi:i:37\tXf:f:" xstr(PI) "\tXd:d:" xstr(E) "\tXZ:Z:" HELLO "\tXH:H:" BEEF "\tXB:B:c,-2,0,+2\tB0:B:i,-2147483648,-1,0,1,2147483647\tB1:B:I,0,1,2147483648,4294967295\tB2:B:s,-32768,-1,0,1,32767\tB3:B:S,0,1,32768,65535\tB4:B:c,-128,-1,0,1,127\tB5:B:C,0,1,127,255\tBf:B:f,-3.14159,2.71828\tZZ:i:1000000\tF2:d:2.46801\tY1:i:-2147483648\tY2:i:-2147483647\tY3:i:-1\tY4:i:0\tY5:i:1\tY6:i:2147483647\tY7:i:2147483648\tY8:i:4294967295\n";

    // Canonical form of the alignment record above, as output by sam_format1()
    static const char r1[] = "r1\t0\tone\t500\t20\t8M\t*\t0\t0\tATGCATGC\tqqqqqqqq\tXi:i:37\tXf:f:3.14159\tXd:d:2.71828\tXZ:Z:" NEW_HELLO "\tXH:H:" BEEF "\tXB:B:c,-2,0,2\tB0:B:i,-2147483648,-1,0,1,2147483647\tB1:B:I,0,1,2147483648,4294967295\tB2:B:s,-32768,-1,0,1,32767\tB3:B:S,0,1,32768,65535\tB4:B:c,-128,-1,0,1,127\tB5:B:C,0,1,127,255\tBf:B:f,-3.14159,2.71828\tZZ:i:1000000\tF2:f:9.8765\tY1:i:-2147483648\tY2:i:-2147483647\tY3:i:-1\tY4:i:0\tY5:i:1\tY6:i:2147483647\tY7:i:2147483648\tY8:i:4294967295\tN0:i:-1234\tN1:i:1234\tN2:i:-2\tN3:i:3\tF1:f:4.5678\tN4:B:S,65535,32768,1,0\tN5:i:4242";

    samFile *in = sam_open(sam, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    uint8_t *p;
    kstring_t ks = { 0, 0, NULL };
    int64_t b0vals[5] = { -2147483648LL,-1,0,1,2147483647LL }; // i
    int64_t b1vals[4] = { 0,1,2147483648LL,4294967295LL };     // I
    int64_t b2vals[5] = { -32768,-1,0,1,32767 };           // s
    int64_t b3vals[4] = { 0,1,32768,65535 };               // S
    int64_t b4vals[5] = { -128,-1,0,1,127 };               // c
    int64_t b5vals[4] = { 0,1,127,255 };                   // C
    // NB: Floats not doubles below!
    // See https://randomascii.wordpress.com/2012/06/26/doubles-are-not-floats-so-dont-compare-them/
    float bfvals[2] = { -3.14159f, 2.71828f };

    int8_t n4v1[] = { -128, -64, -32, -16, -8, -4, -2, -1,
                      0, 1, 2, 4, 8, 16, 32, 64, 127 };
    uint32_t n4v2[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1234, 5678, 1U << 31, 0 };
    int16_t n4v3[] = { -32768, -1, 0, 1, 32767 };
    float n4v4[] = { 0, 1, 2, 10, 20, 30, 1.5, -2.5 };
    uint8_t n4v5[] = { 0, 255 };
    int32_t n4v6[] = { -2147483647 - 1, 10, -1, 0, 1, 2147483647 };
    uint16_t n4v7[] = { 65535, 32768, 1, 0 };

    int32_t ival = -1234;
    uint32_t uval = 1234;
    float f1 = 4.5678;
    float f2 = 9.8765;

    size_t nvals, i;

    if (sam_read1(in, header, aln) >= 0) {
        if ((p = check_bam_aux_get(aln, "XA", 'A')) && bam_aux2A(p) != 'k')
            fail("XA field is '%c', expected 'k'", bam_aux2A(p));

        bam_aux_del(aln,p);
        if (bam_aux_get(aln,"XA"))
            fail("XA field was not deleted");

        if ((p = check_bam_aux_get(aln, "Xi", 'C')) && bam_aux2i(p) != 37)
            fail("Xi field is %"PRId64", expected 37", bam_aux2i(p));

        if ((p = check_bam_aux_get(aln, "Xf", 'f')) && fabs(bam_aux2f(p) - PI) > 1E-6)
            fail("Xf field is %.12f, expected pi", bam_aux2f(p));

        if ((p = check_bam_aux_get(aln, "Xd", 'd')) && fabs(bam_aux2f(p) - E) > 1E-6)
            fail("Xf field is %.12f, expected e", bam_aux2f(p));

        if ((p = check_bam_aux_get(aln, "XZ", 'Z')) && strcmp(bam_aux2Z(p), HELLO) != 0)
            fail("XZ field is \"%s\", expected \"%s\"", bam_aux2Z(p), HELLO);

        bam_aux_update_str(aln,"XZ",strlen(NEW_HELLO)+1,NEW_HELLO);
        if ((p = check_bam_aux_get(aln, "XZ", 'Z')) && strcmp(bam_aux2Z(p), NEW_HELLO) != 0)
            fail("XZ field is \"%s\", expected \"%s\"", bam_aux2Z(p), NEW_HELLO);


        if ((p = check_bam_aux_get(aln, "XH", 'H')) && strcmp(bam_aux2Z(p), BEEF) != 0)
            fail("XH field is \"%s\", expected \"%s\"", bam_aux2Z(p), BEEF);

        if ((p = check_bam_aux_get(aln, "XB", 'B'))
            && ! (memcmp(p, "Bc", 2) == 0
                  && memcmp(p + 2, "\x03\x00\x00\x00\xfe\x00\x02", 7) == 0))
            fail("XB field is %c,..., expected c,-2,0,+2", p[1]);

        check_int_B_array(aln, "B0", NELE(b0vals), b0vals);
        check_int_B_array(aln, "B1", NELE(b1vals), b1vals);
        check_int_B_array(aln, "B2", NELE(b2vals), b2vals);
        check_int_B_array(aln, "B3", NELE(b3vals), b3vals);
        check_int_B_array(aln, "B4", NELE(b4vals), b4vals);
        check_int_B_array(aln, "B5", NELE(b5vals), b5vals);

        nvals = NELE(bfvals);
        if ((p = check_bam_aux_get(aln, "Bf", 'B')) != NULL) {
            if (bam_auxB_len(p) != nvals)
                fail("Wrong length reported for Bf field, got %d, expected %zd\n",
                     bam_auxB_len(p), nvals);

            for (i = 0; i < nvals; i++) {
                if (bam_auxB2f(p, i) != bfvals[i]) {
                    fail("Wrong value from bam_auxB2f for Bf field index %zd, "
                         "got %f expected %f\n",
                         i, bam_auxB2f(p, i), bfvals[i]);
                }
            }
        }

        if ((p = check_bam_aux_get(aln, "ZZ", 'I')) && bam_aux2i(p) != 1000000)
            fail("ZZ field is %"PRId64", expected 1000000", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y1")) && bam_aux2i(p) != -2147483647-1)
            fail("Y1 field is %"PRId64", expected -2^31", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y2")) && bam_aux2i(p) != -2147483647)
            fail("Y2 field is %"PRId64", expected -2^31+1", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y3")) && bam_aux2i(p) != -1)
            fail("Y3 field is %"PRId64", expected -1", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y4")) && bam_aux2i(p) != 0)
            fail("Y4 field is %"PRId64", expected 0", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y5")) && bam_aux2i(p) != 1)
            fail("Y5 field is %"PRId64", expected 1", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y6")) && bam_aux2i(p) != 2147483647)
            fail("Y6 field is %"PRId64", expected 2^31-1", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y7")) && bam_aux2i(p) != 2147483648LL)
            fail("Y7 field is %"PRId64", expected 2^31", bam_aux2i(p));

        if ((p = bam_aux_get(aln, "Y8")) && bam_aux2i(p) != 4294967295LL)
            fail("Y8 field is %"PRId64", expected 2^32-1", bam_aux2i(p));

        // Try appending some new tags
        if (bam_aux_append(aln, "N0", 'i', sizeof(ival), (uint8_t *) &ival) != 0)
            fail("Failed to append N0:i tag");

        if ((p = bam_aux_get(aln, "N0")) && bam_aux2i(p) != ival)
            fail("N0 field is %"PRId64", expected %d", bam_aux2i(p), ival);

        if (bam_aux_append(aln, "N1", 'I', sizeof(uval), (uint8_t *) &uval) != 0)
            fail("failed to append N1:I tag");

        if ((p = bam_aux_get(aln, "N1")) && bam_aux2i(p) != uval)
            fail("N1 field is %"PRId64", expected %u", bam_aux2i(p), uval);

        // Append tags with bam_aux_update_int()
        if (bam_aux_update_int(aln, "N2", -2) < 0)
            fail("failed to append N2:c tag");

        if (bam_aux_update_int(aln, "N3", 3) < 0)
            fail("failed to append N3:C tag");

        p = bam_aux_get(aln, "N2");
        if (!p)
            fail("failed to retrieve N2 tag");
        else if (*p != 'c' || bam_aux2i(p) != -2)
            fail("N2 field is %c:%"PRId64", expected c:-2", *p, bam_aux2i(p));

        p = bam_aux_get(aln, "N3");
        if (!p)
            fail("failed to retrieve N3 tag");
        else if (*p != 'C' || bam_aux2i(p) != 3)
            fail("N3 field is %c:%"PRId64", expected C:3", *p, bam_aux2i(p));

        // Try changing values with bam_aux_update_int()
        i = test_update_int(aln, "N2", 2, 'C', "N3", 3, 'C');
        if (i == 0) test_update_int(aln, "N2", 1234, 'S', "N3", 3, 'C');
        if (i == 0) test_update_int(aln, "N2", -1, 's', "N3", 3, 'C');
        if (i == 0) test_update_int(aln, "N2", 4294967295U, 'I', "N3", 3, 'C');
        if (i == 0) test_update_int(aln, "N2", -2, 'i', "N3", 3, 'C');

        // Append a value with bam_aux_update_float()
        if (bam_aux_update_float(aln, "F1", f1) < 0)
            fail("append F1:f tag");

        p = bam_aux_get(aln, "F1");
        if (!p)
            fail("retrieve F1 tag");
        else if (*p != 'f' || bam_aux2f(p) != f1)
            fail("F1 field is %c:%e, expected f:%e", *p, bam_aux2f(p), f1);

        // Change a double tag to a float
        if (bam_aux_update_float(aln, "F2", f2) < 0)
            fail("update F2 tag");

        p = bam_aux_get(aln, "F2");
        if (!p)
            fail("retrieve F2 tag");
        else if (*p != 'f' || bam_aux2f(p) != f2)
            fail("F2 field is %c:%e, expected f:%e", *p, bam_aux2f(p), f2);

        // Check the next one is intact too
        p = bam_aux_get(aln, "Y1");
        if (!p)
            fail("retrieve Y1 tag");
        else if (*p != 'i' && bam_aux2i(p) != -2147483647-1)
            fail("Y1 field is %"PRId64", expected -2^31", bam_aux2i(p));

        // bam_aux_update_array tests
        // append a new array
        i = test_update_array(aln, "N4", 'c', NELE(n4v1), n4v1, "\0\0", 0, 0);

        // Add a sentinal to check resizes work
        if (i == 0) i = test_update_int(aln, "N5", 4242, 'S', "\0\0", 0, 0);

        // alter the array tag a few times
        if (i == 0)
            i = test_update_array(aln, "N4", 'I', NELE(n4v2), n4v2,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 's', NELE(n4v3), n4v3,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 'f', NELE(n4v4), n4v4,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 'c', NELE(n4v5), n4v5,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 'i', NELE(n4v6), n4v6,
                                  "N5", 4242, 'S');
        if (i == 0)
            i = test_update_array(aln, "N4", 'S', NELE(n4v7), n4v7,
                                  "N5", 4242, 'S');

        if (sam_format1(header, aln, &ks) < 0)
            fail("can't format record");

        if (strcmp(ks.s, r1) != 0)
            fail("record formatted incorrectly: \"%s\"", ks.s);

        free(ks.s);
    }
    else fail("can't read record");

    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);

    return 1;
}

static void iterators1(void)
{
    hts_itr_destroy(sam_itr_queryi(NULL, HTS_IDX_REST, 0, 0));
    hts_itr_destroy(sam_itr_queryi(NULL, HTS_IDX_NONE, 0, 0));
}

static void copy_check_alignment(const char *infname, const char *informat,
    const char *outfname, const char *outmode, const char *outref)
{
    samFile *in = sam_open(infname, "r");
    samFile *out = sam_open(outfname, outmode);
    bam1_t *aln = bam_init1();
    bam_hdr_t *header = NULL;
    int res;

    if (!in) {
        fail("couldn't open %s", infname);
        goto err;
    }
    if (!out) {
        fail("couldn't open %s with mode %s", outfname, outmode);
        goto err;
    }
    if (!aln) {
        fail("bam_init1() failed");
        goto err;
    }

    if (outref) {
        if (hts_set_opt(out, CRAM_OPT_REFERENCE, outref) < 0) {
            fail("setting reference %s for %s", outref, outfname);
            goto err;
        }
    }

    header = sam_hdr_read(in);
    if (!header) {
        fail("reading header from %s", infname);
        goto err;
    }
    if (sam_hdr_write(out, header) < 0) fail("writing headers to %s", outfname);

    while ((res = sam_read1(in, header, aln)) >= 0) {
        int mod4 = ((intptr_t) bam_get_cigar(aln)) % 4;
        if (mod4 != 0)
            fail("%s CIGAR not 4-byte aligned; offset is 4k+%d for \"%s\"",
                 informat, mod4, bam_get_qname(aln));

        if (sam_write1(out, header, aln) < 0) fail("writing to %s", outfname);
    }
    if (res < -1) {
        fail("failed to read alignment from %s", infname);
    }

 err:
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    if (in) sam_close(in);
    if (out) sam_close(out);
}

static void samrecord_layout(void)
{
    static const char qnames[] = "data:,"
"@SQ\tSN:CHROMOSOME_II\tLN:5000\n"
    "a\t0\tCHROMOSOME_II\t100\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
   "bc\t0\tCHROMOSOME_II\t200\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
  "def\t0\tCHROMOSOME_II\t300\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
 "ghij\t0\tCHROMOSOME_II\t400\t10\t4M\t*\t0\t0\tATGC\tqqqq\n"
"klmno\t0\tCHROMOSOME_II\t500\t10\t4M\t*\t0\t0\tATGC\tqqqq\n";

    size_t bam1_t_size, bam1_t_size2;

    bam1_t_size = 36 + sizeof (int) + 4 + sizeof (char *);
#ifndef BAM_NO_ID
    bam1_t_size += 8;
#endif
    bam1_t_size2 = bam1_t_size + 4;  // Account for padding on some platforms

    if (sizeof (bam1_core_t) != 36)
        fail("sizeof bam1_core_t is %zu, expected 36", sizeof (bam1_core_t));

    if (sizeof (bam1_t) != bam1_t_size && sizeof (bam1_t) != bam1_t_size2)
        fail("sizeof bam1_t is %zu, expected either %zu or %zu",
             sizeof(bam1_t), bam1_t_size, bam1_t_size2);

    copy_check_alignment(qnames, "SAM",
                         "test/sam_alignment.tmp.bam", "wb", NULL);
    copy_check_alignment("test/sam_alignment.tmp.bam", "BAM",
                         "test/sam_alignment.tmp.cram", "wc", "test/ce.fa");
    copy_check_alignment("test/sam_alignment.tmp.cram", "CRAM",
                         "test/sam_alignment.tmp.sam_", "w", NULL);
}

static void faidx1(const char *filename)
{
    int n, n_exp = 0, n_fq_exp = 0;
    char tmpfilename[FILENAME_MAX], line[500];
    FILE *fin, *fout;
    faidx_t *fai;

    fin = fopen(filename, "rb");
    if (fin == NULL) fail("can't open %s\n", filename);
    sprintf(tmpfilename, "%s.tmp", filename);
    fout = fopen(tmpfilename, "wb");
    if (fout == NULL) fail("can't create temporary %s\n", tmpfilename);
    while (fgets(line, sizeof line, fin)) {
        if (line[0] == '>') n_exp++;
        if (line[0] == '+' && line[1] == '\n') n_fq_exp++;
        fputs(line, fout);
    }
    fclose(fin);
    fclose(fout);

    if (n_exp == 0 && n_fq_exp != 0) {
        // probably a fastq file
        n_exp = n_fq_exp;
    }

    if (fai_build(tmpfilename) < 0) fail("can't index %s", tmpfilename);
    fai = fai_load(tmpfilename);
    if (fai == NULL) { fail("can't load faidx file %s", tmpfilename); return; }

    n = faidx_fetch_nseq(fai);
    if (n != n_exp)
        fail("%s: faidx_fetch_nseq returned %d, expected %d", filename, n, n_exp);

    n = faidx_nseq(fai);
    if (n != n_exp)
        fail("%s: faidx_nseq returned %d, expected %d", filename, n, n_exp);

    fai_destroy(fai);
}

static void check_enum1(void)
{
    // bgzf_compression() returns int, but enjoys this correspondence
    if (no_compression != 0) fail("no_compression is %d", no_compression);
    if (gzip != 1) fail("gzip is %d", gzip);
    if (bgzf != 2) fail("bgzf is %d", bgzf);
}

int main(int argc, char **argv)
{
    int i;

    status = EXIT_SUCCESS;

    aux_fields1();
    iterators1();
    samrecord_layout();
    check_enum1();
    for (i = 1; i < argc; i++) faidx1(argv[i]);

    return status;
}
