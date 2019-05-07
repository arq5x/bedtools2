/*
Copyright (c) 2012-2013 Genome Research Ltd.
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

#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "cram/cram.h"
#include "cram/os.h"

cram_stats *cram_stats_create(void) {
    return calloc(1, sizeof(cram_stats));
}

void cram_stats_add(cram_stats *st, int32_t val) {
    st->nsamp++;

    //assert(val >= 0);

    if (val < MAX_STAT_VAL && val >= 0) {
        st->freqs[val]++;
    } else {
        khint_t k;
        int r;

        if (!st->h) {
            st->h = kh_init(m_i2i);
        }

        k = kh_put(m_i2i, st->h, val, &r);
        if (r == 0)
            kh_val(st->h, k)++;
        else if (r != -1)
            kh_val(st->h, k) = 1;
        else
            ; // FIXME: handle error
    }
}

void cram_stats_del(cram_stats *st, int32_t val) {
    st->nsamp--;

    //assert(val >= 0);

    if (val < MAX_STAT_VAL && val >= 0) {
        st->freqs[val]--;
        assert(st->freqs[val] >= 0);
    } else if (st->h) {
        khint_t k = kh_get(m_i2i, st->h, val);

        if (k != kh_end(st->h)) {
            if (--kh_val(st->h, k) == 0)
                kh_del(m_i2i, st->h, k);
        } else {
            hts_log_warning("Failed to remove val %d from cram_stats", val);
            st->nsamp++;
        }
    } else {
        hts_log_warning("Failed to remove val %d from cram_stats", val);
        st->nsamp++;
    }
}

#if DEBUG_CRAM_STATS
void cram_stats_dump(cram_stats *st) {
    int i;
    fprintf(stderr, "cram_stats:\n");
    for (i = 0; i < MAX_STAT_VAL; i++) {
        if (!st->freqs[i])
            continue;
        fprintf(stderr, "\t%d\t%d\n", i, st->freqs[i]);
    }
    if (st->h) {
        khint_t k;
        for (k = kh_begin(st->h); k != kh_end(st->h); k++) {
            if (!kh_exist(st->h, k))
                continue;

            fprintf(stderr, "\t%d\t%d\n", kh_key(st->h, k), kh_val(st->h, k));
        }
    }
}
#endif

/*
 * Computes entropy from integer frequencies for various encoding methods and
 * picks the best encoding.
 *
 * FIXME: we could reuse some of the code here for the actual encoding
 * parameters too. Eg the best 'k' for SUBEXP or the code lengths for huffman.
 *
 * Returns the best codec to use.
 */
enum cram_encoding cram_stats_encoding(cram_fd *fd, cram_stats *st) {
    int nvals, i, ntot = 0, max_val = 0, min_val = INT_MAX;
    int *vals = NULL, *freqs = NULL, vals_alloc = 0;

#if DEBUG_CRAM_STATS
    cram_stats_dump(st);
#endif

    /* Count number of unique symbols */
    for (nvals = i = 0; i < MAX_STAT_VAL; i++) {
        if (!st->freqs[i])
            continue;
        if (nvals >= vals_alloc) {
            vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
            vals  = realloc(vals,  vals_alloc * sizeof(int));
            freqs = realloc(freqs, vals_alloc * sizeof(int));
            if (!vals || !freqs) {
                if (vals)  free(vals);
                if (freqs) free(freqs);
                return E_HUFFMAN; // Cannot do much else atm
            }
        }
        vals[nvals] = i;
        freqs[nvals] = st->freqs[i];
        ntot += freqs[nvals];
        if (max_val < i) max_val = i;
        if (min_val > i) min_val = i;
        nvals++;
    }
    if (st->h) {
        khint_t k;
        int i;

        for (k = kh_begin(st->h); k != kh_end(st->h); k++) {
            if (!kh_exist(st->h, k))
                continue;

            if (nvals >= vals_alloc) {
                vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
                vals  = realloc(vals,  vals_alloc * sizeof(int));
                freqs = realloc(freqs, vals_alloc * sizeof(int));
                if (!vals || !freqs)
                    return E_HUFFMAN; // Cannot do much else atm
            }
            i = kh_key(st->h, k);
            vals[nvals]=i;
            freqs[nvals] = kh_val(st->h, k);
            ntot += freqs[nvals];
            if (max_val < i) max_val = i;
            if (min_val > i) min_val = i;
            nvals++;
        }
    }

    st->nvals = nvals;
    assert(ntot == st->nsamp);

    free(vals);
    free(freqs);

    /*
     * Simple policy that everything is external unless it can be
     * encoded using zero bits as a unary item huffman table.
     */
    return nvals <= 1 ? E_HUFFMAN : E_EXTERNAL;
}

void cram_stats_free(cram_stats *st) {
    if (st->h)
        kh_destroy(m_i2i, st->h);
    free(st);
}
