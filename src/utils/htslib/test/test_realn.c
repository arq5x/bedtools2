/*  test/test_realn.c -- test sam_prob_realn() function

    Copyright (C) 2018 Genome Research Ltd.

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

#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"

void usage(const char *prog) {
    fprintf(stderr, "Usage: %s -i <in.sam> -o <out.sam> -f <ref.fa>\n", prog);
}

int main(int argc, char **argv) {
    htsFile *in = NULL;
    htsFile *out = NULL;
    char *in_name = "-";
    char *out_name = "-";
    char *ref_name = NULL;
    char *ref_seq = NULL;
    char modew[8] = "w";
    faidx_t *fai = NULL;
    bam_hdr_t *hdr = NULL;
    bam1_t *rec = NULL;
    int c, res, last_ref = -1, ref_len = 0;
    int adjust = 0, extended = 0, recalc = 0, flags = 0;

    while ((c = getopt(argc, argv, "aef:hi:o:r")) >= 0) {
        switch (c) {
        case 'a': adjust = 1; break;
        case 'e': extended = 1; break;
        case 'f': ref_name = optarg; break;
        case 'h': usage(argv[0]); return EXIT_SUCCESS;
        case 'i': in_name = optarg; break;
        case 'o': out_name = optarg; break;
        case 'r': recalc = 1; break;
        default: usage(argv[0]); return EXIT_FAILURE;
        }
    }

    if (!ref_name) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    flags = (adjust ? 1 : 0) | (extended ? 2 : 0) | (recalc ? 4 : 0);

    fai = fai_load(ref_name);
    if (!fai) {
        fprintf(stderr, "Couldn't load reference %s\n", ref_name);
        goto fail;
    }

    rec = bam_init1();
    if (!rec) {
        perror(NULL);
        goto fail;
    }

    in = hts_open(in_name, "r");
    if (!in) {
        fprintf(stderr, "Couldn't open %s : %s\n", in_name, strerror(errno));
        goto fail;
    }

    hdr = sam_hdr_read(in);
    if (!hdr) {
        fprintf(stderr, "Couldn't read header for %s\n", in_name);
        goto fail;
    }

    out = hts_open(out_name, modew);
    if (!out) {
        fprintf(stderr, "Couldn't open %s : %s\n", out_name, strerror(errno));
        goto fail;
    }

    if (sam_hdr_write(out, hdr) < 0) {
        fprintf(stderr, "Couldn't write header to %s : %s\n",
                out_name, strerror(errno));
        goto fail;
    }

    while ((res = sam_read1(in, hdr, rec)) >= 0) {
        if (rec->core.tid >= hdr->n_targets) {
            fprintf(stderr, "Invalid BAM reference id %d\n", rec->core.tid);
            goto fail;
        }
        if (last_ref != rec->core.tid && rec->core.tid >= 0) {
            free(ref_seq);
            ref_seq = faidx_fetch_seq(fai, hdr->target_name[rec->core.tid],
                                      0, INT_MAX, &ref_len);
            if (!ref_seq) {
                fprintf(stderr, "Couldn't get reference %s\n",
                        hdr->target_name[rec->core.tid]);
                goto fail;
            }
            last_ref = rec->core.tid;
        }
        if (rec->core.tid >= 0) {
            res = sam_prob_realn(rec, ref_seq, ref_len, flags);
            if (res <= -4) {
                fprintf(stderr, "Error running sam_prob_realn : %s\n",
                        strerror(errno));
                goto fail;
            }
        }
        if (sam_write1(out, hdr, rec) < 0) {
            fprintf(stderr, "Error writing to %s\n", out_name);
            goto fail;
        }
    }
    res = hts_close(in);
    in = NULL;
    if (res < 0) {
        fprintf(stderr, "Error closing %s\n", in_name);
        goto fail;
    }

    res = hts_close(out);
    out = NULL;
    if (res < 0) {
        fprintf(stderr, "Error closing %s\n", out_name);
        goto fail;
    }

    bam_hdr_destroy(hdr);
    bam_destroy1(rec);
    free(ref_seq);
    fai_destroy(fai);

    return EXIT_SUCCESS;

 fail:
    if (hdr) bam_hdr_destroy(hdr);
    if (rec) bam_destroy1(rec);
    if (in) hts_close(in);
    if (out) hts_close(out);
    free(ref_seq);
    fai_destroy(fai);
    return EXIT_FAILURE;
}
