/*  faidx.c -- FASTA random access.

    Copyright (C) 2008, 2009, 2013-2017 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.

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

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <assert.h>

#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/hfile.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "hts_internal.h"

typedef struct {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

struct __faidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
};

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline int fai_insert_index(faidx_t *idx, const char *name, int64_t len, int line_len, int line_blen, uint64_t offset)
{
    if (!name) {
        hts_log_error("Malformed line");
        return -1;
    }

    char *name_key = strdup(name);
    int absent;
    khint_t k = kh_put(s, idx->hash, name_key, &absent);
    faidx1_t *v = &kh_value(idx->hash, k);

    if (! absent) {
        hts_log_warning("Ignoring duplicate sequence \"%s\" at byte offset %"PRIu64"", name, offset);
        free(name_key);
        return 0;
    }

    if (idx->n == idx->m) {
        char **tmp;
        idx->m = idx->m? idx->m<<1 : 16;
        if (!(tmp = (char**)realloc(idx->name, sizeof(char*) * idx->m))) {
            hts_log_error("Out of memory");
            return -1;
        }
        idx->name = tmp;
    }
    idx->name[idx->n++] = name_key;
    v->len = len;
    v->line_len = line_len;
    v->line_blen = line_blen;
    v->offset = offset;

    return 0;
}

faidx_t *fai_build_core(BGZF *bgzf)
{
    kstring_t name = { 0, 0, NULL };
    int c;
    int line_len, line_blen, state;
    int l1, l2;
    faidx_t *idx;
    uint64_t offset;
    int64_t len;

    idx = (faidx_t*)calloc(1, sizeof(faidx_t));
    idx->hash = kh_init(s);
    len = line_len = line_blen = -1; state = 0; l1 = l2 = -1; offset = 0;
    while ( (c=bgzf_getc(bgzf))>=0 ) {
        if (c == '\n') { // an empty line
            if (state == 1) {
                offset = bgzf_utell(bgzf);
                continue;
            } else if ((state == 0 && len < 0) || state == 2) continue;
            else if (state == 0) { state = 2; continue; }
        }
        if (c == '>') { // fasta header
            if (len >= 0) {
                if (fai_insert_index(idx, name.s, len, line_len, line_blen, offset) != 0)
                    goto fail;
            }

            name.l = 0;
            while ((c = bgzf_getc(bgzf)) >= 0)
                if (! isspace(c)) kputc_(c, &name);
                else if (name.l > 0 || c == '\n') break;
            kputsn("", 0, &name);

            if ( c<0 ) {
                hts_log_error("The last entry has no sequence");
                goto fail;
            }
            if (c != '\n') while ( (c=bgzf_getc(bgzf))>=0 && c != '\n');
            state = 1; len = 0;
            offset = bgzf_utell(bgzf);
        } else {
            if (state == 3) {
                hts_log_error("Inlined empty line is not allowed in sequence '%s'", name.s);
                goto fail;
            }
            if (state == 2) state = 3;
            l1 = l2 = 0;
            do {
                ++l1;
                if (isgraph(c)) ++l2;
            } while ( (c=bgzf_getc(bgzf))>=0 && c != '\n');
            if (state == 3 && l2) {
                hts_log_error("Different line length in sequence '%s'", name.s);
                goto fail;
            }
            ++l1; len += l2;
            if (state == 1) line_len = l1, line_blen = l2, state = 0;
            else if (state == 0) {
                if (l1 != line_len || l2 != line_blen) state = 2;
            }
        }
    }

    if (len >= 0) {
        if (fai_insert_index(idx, name.s, len, line_len, line_blen, offset) != 0)
            goto fail;
    } else {
        goto fail;
    }

    free(name.s);
    return idx;

fail:
    free(name.s);
    fai_destroy(idx);
    return NULL;
}

static int fai_save(const faidx_t *fai, hFILE *fp) {
    khint_t k;
    int i;
    char buf[96]; // Must be big enough for format below.

    for (i = 0; i < fai->n; ++i) {
        faidx1_t x;
        k = kh_get(s, fai->hash, fai->name[i]);
        assert(k < kh_end(fai->hash));
        x = kh_value(fai->hash, k);
        snprintf(buf, sizeof(buf),
                 "\t%"PRId64"\t%"PRIu64"\t%"PRId32"\t%"PRId32"\n",
                 x.len, x.offset, x.line_blen, x.line_len);
        if (hputs(fai->name[i], fp) != 0) return -1;
        if (hputs(buf, fp) != 0) return -1;
    }
    return 0;
}

static faidx_t *fai_read(hFILE *fp, const char *fname)
{
    faidx_t *fai;
    char *buf = NULL, *p;
    int line_len, line_blen, n;
    int64_t len;
    uint64_t offset;
    ssize_t l, lnum = 1;

    fai = (faidx_t*)calloc(1, sizeof(faidx_t));
    if (!fai) return NULL;

    fai->hash = kh_init(s);
    if (!fai->hash) goto fail;

    buf = (char*)calloc(0x10000, 1);
    if (!buf) goto fail;

    while ((l = hgetln(buf, 0x10000, fp)) > 0) {
        for (p = buf; *p && !isspace_c(*p); ++p);
        if (p - buf < l) {
            *p = 0; ++p;
        }
        n = sscanf(p, "%"SCNd64"%"SCNu64"%d%d", &len, &offset, &line_blen, &line_len);
        if (n != 4) {
            hts_log_error("Could not understand FAI %s line %zd", fname, lnum);
            goto fail;
        }
        if (fai_insert_index(fai, buf, len, line_len, line_blen, offset) != 0) {
            goto fail;
        }
        if (buf[l - 1] == '\n') ++lnum;
    }

    if (l < 0) {
        hts_log_error("Error while reading %s: %s", fname, strerror(errno));
        goto fail;
    }
    free(buf);
    return fai;

 fail:
    free(buf);
    fai_destroy(fai);
    return NULL;
}

void fai_destroy(faidx_t *fai)
{
    int i;
    if (!fai) return;
    for (i = 0; i < fai->n; ++i) free(fai->name[i]);
    free(fai->name);
    kh_destroy(s, fai->hash);
    if (fai->bgzf) bgzf_close(fai->bgzf);
    free(fai);
}

int fai_build3(const char *fn, const char *fnfai, const char *fngzi)
{
    kstring_t fai_kstr = { 0, 0, NULL };
    kstring_t gzi_kstr = { 0, 0, NULL };
    BGZF *bgzf = NULL;
    hFILE *fp = NULL;
    faidx_t *fai = NULL;
    int save_errno, res;

    if (!fnfai) {
        if (ksprintf(&fai_kstr, "%s.fai", fn) < 0) goto fail;
        fnfai = fai_kstr.s;
    }
    if (!fngzi) {
        if (ksprintf(&gzi_kstr, "%s.gzi", fn) < 0) goto fail;
        fngzi = gzi_kstr.s;
    }

    bgzf = bgzf_open(fn, "r");
    if ( !bgzf ) {
        hts_log_error("Failed to open the FASTA file %s", fn);
        goto fail;
    }
    if ( bgzf->is_compressed ) {
        if (bgzf_index_build_init(bgzf) != 0) {
            hts_log_error("Failed to allocate bgzf index");
            goto fail;
        }
    }
    fai = fai_build_core(bgzf);
    if ( !fai ) {
        if (bgzf->is_compressed && bgzf->is_gzip) {
            hts_log_error("Cannot index files compressed with gzip, please use bgzip");
        }
        goto fail;
    }
    if ( bgzf->is_compressed ) {
        if (bgzf_index_dump(bgzf, fngzi, NULL) < 0) {
            hts_log_error("Failed to make bgzf index %s", fngzi);
            goto fail;
        }
    }
    res = bgzf_close(bgzf);
    bgzf = NULL;
    if (res < 0) {
        hts_log_error("Error on closing %s : %s", fn, strerror(errno));
        goto fail;
    }
    fp = hopen(fnfai, "wb");
    if ( !fp ) {
        hts_log_error("Failed to open FASTA index %s : %s", fnfai, strerror(errno));
        goto fail;
    }
    if (fai_save(fai, fp) != 0) {
        hts_log_error("Failed to write FASTA index %s : %s", fnfai, strerror(errno));
        goto fail;
    }
    if (hclose(fp) != 0) {
        hts_log_error("Failed on closing FASTA index %s : %s", fnfai, strerror(errno));
        goto fail;
    }

    free(fai_kstr.s);
    free(gzi_kstr.s);
    fai_destroy(fai);
    return 0;

 fail:
    save_errno = errno;
    free(fai_kstr.s);
    free(gzi_kstr.s);
    bgzf_close(bgzf);
    fai_destroy(fai);
    errno = save_errno;
    return -1;
}

int fai_build(const char *fn) {
    return fai_build3(fn, NULL, NULL);
}

faidx_t *fai_load3(const char *fn, const char *fnfai, const char *fngzi,
                   int flags)
{
    kstring_t fai_kstr = { 0, 0, NULL };
    kstring_t gzi_kstr = { 0, 0, NULL };
    hFILE *fp = NULL;
    faidx_t *fai = NULL;
    int res;

    if (fn == NULL)
        return NULL;

    if (fnfai == NULL) {
        if (ksprintf(&fai_kstr, "%s.fai", fn) < 0) goto fail;
        fnfai = fai_kstr.s;
    }
    if (fngzi == NULL) {
        if (ksprintf(&gzi_kstr, "%s.gzi", fn) < 0) goto fail;
        fngzi = gzi_kstr.s;
    }

    fp = hopen(fnfai, "rb");

    if (fp == 0) {
        if (!(flags & FAI_CREATE) || errno != ENOENT) {
            hts_log_error("Failed to open FASTA index %s: %s", fnfai, strerror(errno));
            goto fail;
        }

        hts_log_info("Build FASTA index");

        if (fai_build3(fn, fnfai, fngzi) < 0) {
            goto fail;
        }

        fp = hopen(fnfai, "rb");
        if (fp == 0) {
            hts_log_error("Failed to open FASTA index %s: %s", fnfai, strerror(errno));
            goto fail;
        }
    }

    fai = fai_read(fp, fnfai);
    if (fai == NULL) {
        hts_log_error("Failed to read FASTA index %s", fnfai);
        goto fail;
    }

    res = hclose(fp);
    fp = NULL;
    if (res < 0) {
        hts_log_error("Failed on closing FASTA index %s : %s", fnfai, strerror(errno));
        goto fail;
    }

    fai->bgzf = bgzf_open(fn, "rb");
    if (fai->bgzf == 0) {
        hts_log_error("Failed to open FASTA file %s", fn);
        goto fail;
    }
    if ( fai->bgzf->is_compressed==1 ) {
        if ( bgzf_index_load(fai->bgzf, fngzi, NULL) < 0 ) {
            hts_log_error("Failed to load .gzi index: %s", fngzi);
            goto fail;
        }
    }
    free(fai_kstr.s);
    free(gzi_kstr.s);
    return fai;

 fail:
    if (fai) fai_destroy(fai);
    if (fp) hclose_abruptly(fp);
    free(fai_kstr.s);
    free(gzi_kstr.s);
    return NULL;
}

faidx_t *fai_load(const char *fn)
{
    return fai_load3(fn, NULL, NULL, FAI_CREATE);
}

static char *fai_retrieve(const faidx_t *fai, const faidx1_t *val,
                          long beg, long end, int *len) {
    char *s;
    size_t l;
    int c = 0;
    int ret = bgzf_useek(fai->bgzf,
                         val->offset
                         + beg / val->line_blen * val->line_len
                         + beg % val->line_blen, SEEK_SET);

    if (ret < 0) {
        *len = -1;
        hts_log_error("Failed to retrieve block. (Seeking in a compressed, .gzi unindexed, file?)");
        return NULL;
    }

    l = 0;
    s = (char*)malloc((size_t) end - beg + 2);
    if (!s) {
        *len = -1;
        return NULL;
    }

    while ( l < end - beg && (c=bgzf_getc(fai->bgzf))>=0 )
        if (isgraph(c)) s[l++] = c;
    if (c < 0) {
        hts_log_error("Failed to retrieve block: %s",
            c == -1 ? "unexpected end of file" : "error reading file");
        free(s);
        *len = -1;
        return NULL;
    }

    s[l] = '\0';
    *len = l < INT_MAX ? l : INT_MAX;
    return s;
}

char *fai_fetch(const faidx_t *fai, const char *str, int *len)
{
    char *s, *ep;
    size_t i, l, k, name_end;
    khiter_t iter;
    faidx1_t val;
    khash_t(s) *h;
    long beg, end;

    beg = end = -1;
    h = fai->hash;
    name_end = l = strlen(str);
    s = (char*)malloc(l+1);
    if (!s) {
        *len = -1;
        return NULL;
    }

    // remove space
    for (i = k = 0; i < l; ++i)
        if (!isspace_c(str[i])) s[k++] = str[i];
    s[k] = 0;
    name_end = l = k;
    // determine the sequence name
    for (i = l; i > 0; --i) if (s[i - 1] == ':') break; // look for colon from the end
    if (i > 0) name_end = i - 1;
    if (name_end < l) { // check if this is really the end
        int n_hyphen = 0;
        for (i = name_end + 1; i < l; ++i) {
            if (s[i] == '-') ++n_hyphen;
            else if (!isdigit_c(s[i]) && s[i] != ',') break;
        }
        if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
        s[name_end] = 0;
        iter = kh_get(s, h, s);
        if (iter == kh_end(h)) { // cannot find the sequence name
            iter = kh_get(s, h, str); // try str as the name
            if (iter != kh_end(h)) {
                s[name_end] = ':';
                name_end = l;
            }
        }
    } else iter = kh_get(s, h, str);
    if(iter == kh_end(h)) {
        hts_log_warning("Reference %s not found in FASTA file, returning empty sequence", str);
        free(s);
        *len = -2;
        return 0;
    }
    val = kh_value(h, iter);
    // parse the interval
    if (name_end < l) {
        int save_errno = errno;
        errno = 0;
        for (i = k = name_end + 1; i < l; ++i)
            if (s[i] != ',') s[k++] = s[i];
        s[k] = 0;
        if (s[name_end + 1] == '-') {
            beg = 0;
            i = name_end + 2;
        } else {
            beg = strtol(s + name_end + 1, &ep, 10);
            for (i = ep - s; i < k;) if (s[i++] == '-') break;
        }
        end = i < k? strtol(s + i, &ep, 10) : val.len;
        if (beg > 0) --beg;
        // Check for out of range numbers.  Only going to be a problem on
        // 32-bit platforms with >2Gb sequence length.
        if (errno == ERANGE && (uint64_t) val.len > LONG_MAX) {
            hts_log_error("Positions in range %s are too large for this platform", s);
            free(s);
            *len = -2;
            return NULL;
        }
        errno = save_errno;
    } else beg = 0, end = val.len;
    if (beg >= val.len) beg = val.len;
    if (end >= val.len) end = val.len;
    if (beg > end) beg = end;
    free(s);

    // now retrieve the sequence
    return fai_retrieve(fai, &val, beg, end, len);
}

int faidx_fetch_nseq(const faidx_t *fai)
{
    return fai->n;
}

int faidx_nseq(const faidx_t *fai)
{
    return fai->n;
}

const char *faidx_iseq(const faidx_t *fai, int i)
{
    return fai->name[i];
}

int faidx_seq_len(const faidx_t *fai, const char *seq)
{
    khint_t k = kh_get(s, fai->hash, seq);
    if ( k == kh_end(fai->hash) ) return -1;
    return kh_val(fai->hash, k).len;
}

char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
    khiter_t iter;
    faidx1_t val;

    // Adjust position
    iter = kh_get(s, fai->hash, c_name);
    if (iter == kh_end(fai->hash))
    {
        *len = -2;
        hts_log_error("The sequence \"%s\" not found", c_name);
        return NULL;
    }
    val = kh_value(fai->hash, iter);
    if(p_end_i < p_beg_i) p_beg_i = p_end_i;
    if(p_beg_i < 0) p_beg_i = 0;
    else if(val.len <= p_beg_i) p_beg_i = val.len - 1;
    if(p_end_i < 0) p_end_i = 0;
    else if(val.len <= p_end_i) p_end_i = val.len - 1;

    // Now retrieve the sequence
    return fai_retrieve(fai, &val, p_beg_i, (long) p_end_i + 1, len);
}

int faidx_has_seq(const faidx_t *fai, const char *seq)
{
    khiter_t iter = kh_get(s, fai->hash, seq);
    if (iter == kh_end(fai->hash)) return 0;
    return 1;
}

