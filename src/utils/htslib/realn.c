/*  realn.c -- BAQ calculation and realignment.

    Copyright (C) 2009-2011, 2014-2015 Genome Research Ltd.
    Portions copyright (C) 2009-2011 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <errno.h>
#include <assert.h>
#include "htslib/hts.h"
#include "htslib/sam.h"

int sam_cap_mapq(bam1_t *b, const char *ref, int ref_len, int thres)
{
    uint8_t *seq = bam_get_seq(b), *qual = bam_get_qual(b);
    uint32_t *cigar = bam_get_cigar(b);
    bam1_core_t *c = &b->core;
    int i, x, y, mm, q, len, clip_l, clip_q;
    double t;
    if (thres < 0) thres = 40; // set the default
    mm = q = len = clip_l = clip_q = 0;
    for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
        int j, l = cigar[i]>>4, op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (j = 0; j < l; ++j) {
                int c1, c2, z = y + j;
                if (x+j >= ref_len || ref[x+j] == '\0') break; // out of bounds
                c1 = bam_seqi(seq, z), c2 = seq_nt16_table[(int)ref[x+j]];
                if (c2 != 15 && c1 != 15 && qual[z] >= 13) { // not ambiguous
                    ++len;
                    if (c1 && c1 != c2 && qual[z] >= 13) { // mismatch
                        ++mm;
                        q += qual[z] > 33? 33 : qual[z];
                    }
                }
            }
            if (j < l) break;
            x += l; y += l; len += l;
        } else if (op == BAM_CDEL) {
            for (j = 0; j < l; ++j)
                if (x+j >= ref_len || ref[x+j] == '\0') break;
            if (j < l) break;
            x += l;
        } else if (op == BAM_CSOFT_CLIP) {
            for (j = 0; j < l; ++j) clip_q += qual[y+j];
            clip_l += l;
            y += l;
        } else if (op == BAM_CHARD_CLIP) {
            clip_q += 13 * l;
            clip_l += l;
        } else if (op == BAM_CINS) y += l;
        else if (op == BAM_CREF_SKIP) x += l;
    }
    for (i = 0, t = 1; i < mm; ++i)
        t *= (double)len / (i+1);
    t = q - 4.343 * log(t) + clip_q / 5.;
    if (t > thres) return -1;
    if (t < 0) t = 0;
    t = sqrt((thres - t) / thres) * thres;
    //fprintf(stderr, "%s %lf %d\n", bam_get_qname(b), t, q);
    return (int)(t + .499);
}

static int realn_check_tag(const uint8_t *tg, enum htsLogLevel severity,
                           const char *type, const bam1_t *b) {
    if (*tg != 'Z') {
        hts_log(severity, "Incorrect %s tag type (%c) for read %s",
                type, *tg, bam_get_qname(b));
        return -1;
    }
    if (b->core.l_qseq != strlen((const char *) tg + 1)) {
        hts_log(severity, "Read %s %s tag is wrong length",
                bam_get_qname(b), type);
        return -1;
    }
    return 0;
}

int sam_prob_realn(bam1_t *b, const char *ref, int ref_len, int flag)
{
    int k, i, bw, x, y, yb, ye, xb, xe, apply_baq = flag&1, extend_baq = flag>>1&1, redo_baq = flag&4, fix_bq = 0;
    uint32_t *cigar = bam_get_cigar(b);
    bam1_core_t *c = &b->core;
    probaln_par_t conf = { 0.001, 0.1, 10 };
    uint8_t *bq = NULL, *zq = NULL, *qual = bam_get_qual(b);
    int *state = NULL;
    if ((c->flag & BAM_FUNMAP) || b->core.l_qseq == 0 || qual[0] == (uint8_t)-1)
        return -1; // do nothing

    // test if BQ or ZQ is present, and make sanity checks
    if ((bq = bam_aux_get(b, "BQ")) != NULL) {
        if (!redo_baq) {
            if (realn_check_tag(bq, HTS_LOG_WARNING, "BQ", b) < 0)
                fix_bq = 1;
        }
        ++bq;
    }
    if ((zq = bam_aux_get(b, "ZQ")) != NULL) {
        if (realn_check_tag(zq, HTS_LOG_ERROR, "ZQ", b) < 0)
            return -4;
        ++zq;
    }
    if (bq && redo_baq)
    {
        bam_aux_del(b, bq-1);
        bq = 0;
    }
    if (bq && zq) { // remove the ZQ tag
        bam_aux_del(b, zq-1);
        zq = 0;
    }
    if (!zq && fix_bq) { // Need to fix invalid BQ tag (by realigning)
        assert(bq != NULL);
        bam_aux_del(b, bq-1);
        bq = 0;
    }

    if (bq || zq) {
        if ((apply_baq && zq) || (!apply_baq && bq)) return -3; // in both cases, do nothing
        if (bq && apply_baq) { // then convert BQ to ZQ
            for (i = 0; i < c->l_qseq; ++i)
                qual[i] = qual[i] + 64 < bq[i]? 0 : qual[i] - ((int)bq[i] - 64);
            *(bq - 3) = 'Z';
        } else if (zq && !apply_baq) { // then convert ZQ to BQ
            for (i = 0; i < c->l_qseq; ++i)
                qual[i] += (int)zq[i] - 64;
            *(zq - 3) = 'B';
        }
        return 0;
    }
    // find the start and end of the alignment
    x = c->pos, y = 0, yb = ye = xb = xe = -1;
    for (k = 0; k < c->n_cigar; ++k) {
        int op, l;
        op = cigar[k]&0xf; l = cigar[k]>>4;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (yb < 0) yb = y;
            if (xb < 0) xb = x;
            ye = y + l; xe = x + l;
            x += l; y += l;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
        else if (op == BAM_CDEL) x += l;
        else if (op == BAM_CREF_SKIP) return -1; // do nothing if there is a reference skip
    }
    if (xb == -1) // No matches in CIGAR.
        return -1;
    // set bandwidth and the start and the end
    bw = 7;
    if (abs((xe - xb) - (ye - yb)) > bw)
        bw = abs((xe - xb) - (ye - yb)) + 3;
    conf.bw = bw;
    xb -= yb + bw/2; if (xb < 0) xb = 0;
    xe += c->l_qseq - ye + bw/2;
    if (xe - xb - c->l_qseq > bw)
        xb += (xe - xb - c->l_qseq - bw) / 2, xe -= (xe - xb - c->l_qseq - bw) / 2;
    { // glocal
        uint8_t *seq = bam_get_seq(b);
        uint8_t *tseq; // translated seq A=>0,C=>1,G=>2,T=>3,other=>4
        uint8_t *tref; // translated ref
        uint8_t *q; // Probability of incorrect alignment from probaln_glocal()
        size_t lref = xe > xb ? xe - xb : 1;
        size_t align_lqseq;
        if (extend_baq && lref < c->l_qseq)
            lref = c->l_qseq; // So we can recycle tseq,tref for left,rght below
        // Try to make q,tref,tseq reasonably well aligned
        align_lqseq = ((c->l_qseq + 1) | 0xf) + 1;
        // Overflow check - 3 for *bq, sizeof(int) for *state
        if ((SIZE_MAX - lref) / (3 + sizeof(int)) < align_lqseq) {
            errno = ENOMEM;
            goto fail;
        }

        assert(bq == NULL); // bq was used above, but should now be NULL
        bq = malloc(align_lqseq * 3 + lref);
        if (!bq) goto fail;
        q = bq + align_lqseq;
        tseq = q + align_lqseq;
        tref = tseq + align_lqseq;

        memcpy(bq, qual, c->l_qseq); bq[c->l_qseq] = 0;
        for (i = 0; i < c->l_qseq; ++i)
            tseq[i] = seq_nt16_int[bam_seqi(seq, i)];
        for (i = xb; i < xe; ++i) {
            if (i >= ref_len || ref[i] == '\0') { xe = i; break; }
            tref[i-xb] = seq_nt16_int[seq_nt16_table[(unsigned char)ref[i]]];
        }

        state = malloc(c->l_qseq * sizeof(int));
        if (!state) goto fail;
        if (probaln_glocal(tref, xe-xb, tseq, c->l_qseq, qual,
                           &conf, state, q) == INT_MIN) {
            goto fail;
        }

        if (!extend_baq) { // in this block, bq[] is capped by base quality qual[]
            for (k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
                int op = cigar[k]&0xf, l = cigar[k]>>4;
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    // Sanity check running off the end of the sequence
                    // Can only happen if the alignment is broken
                    if (l > c->l_qseq - y)
                        l = c->l_qseq - y;
                    for (i = y; i < y + l; ++i) {
                        if ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y)) bq[i] = 0;
                        else bq[i] = bq[i] < q[i]? bq[i] : q[i];
                    }
                    x += l; y += l;
                } else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) {
                    // Need sanity check here too.
                    if (l > c->l_qseq - y)
                        l = c->l_qseq - y;
                    y += l;
                } else if (op == BAM_CDEL) {
                    x += l;
                }
            }
            for (i = 0; i < c->l_qseq; ++i) bq[i] = qual[i] - bq[i] + 64; // finalize BQ
        } else { // in this block, bq[] is BAQ that can be larger than qual[] (different from the above!)
            // tseq,tref are no longer needed, so we can steal them to avoid mallocs
            uint8_t *left = tseq;
            uint8_t *rght = tref;
            for (k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
                int op = cigar[k]&0xf, l = cigar[k]>>4;
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    // Sanity check running off the end of the sequence
                    // Can only happen if the alignment is broken
                    if (l > c->l_qseq - y)
                        l = c->l_qseq - y;
                    for (i = y; i < y + l; ++i)
                        bq[i] = ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y))? 0 : q[i];
                    for (left[y] = bq[y], i = y + 1; i < y + l; ++i)
                        left[i] = bq[i] > left[i-1]? bq[i] : left[i-1];
                    for (rght[y+l-1] = bq[y+l-1], i = y + l - 2; i >= y; --i)
                        rght[i] = bq[i] > rght[i+1]? bq[i] : rght[i+1];
                    for (i = y; i < y + l; ++i)
                        bq[i] = left[i] < rght[i]? left[i] : rght[i];
                    x += l; y += l;
                } else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) {
                    // Need sanity check here too.
                    if (l > c->l_qseq - y)
                        l = c->l_qseq - y;
                    y += l;
                } else if (op == BAM_CDEL) {
                    x += l;
                }
            }
            for (i = 0; i < c->l_qseq; ++i) bq[i] = 64 + (qual[i] <= bq[i]? 0 : qual[i] - bq[i]); // finalize BQ
        }
        if (apply_baq) {
            for (i = 0; i < c->l_qseq; ++i) qual[i] -= bq[i] - 64; // modify qual
            bam_aux_append(b, "ZQ", 'Z', c->l_qseq + 1, bq);
        } else bam_aux_append(b, "BQ", 'Z', c->l_qseq + 1, bq);
        free(bq); free(state);
    }
    return 0;

 fail:
    free(bq); free(state);
    return -4;
}
