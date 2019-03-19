/*  faidx.c -- FASTA and FASTQ random access.

    Copyright (C) 2008, 2009, 2013-2018 Genome Research Ltd.
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
    uint32_t line_len, line_blen;
    uint64_t len;
    uint64_t seq_offset;
    uint64_t qual_offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

struct __faidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
    enum fai_format_options format;
};

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline int fai_insert_index(faidx_t *idx, const char *name, uint64_t len, uint32_t line_len, uint32_t line_blen, uint64_t seq_offset, uint64_t qual_offset)
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
        hts_log_warning("Ignoring duplicate sequence \"%s\" at byte offset %"PRIu64"", name, seq_offset);
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
    v->seq_offset = seq_offset;
    v->qual_offset = qual_offset;

    return 0;
}


static faidx_t *fai_build_core(BGZF *bgzf) {
    kstring_t name = { 0, 0, NULL };
    int c, read_done, line_num;
    faidx_t *idx;
    uint64_t seq_offset, qual_offset;
    uint64_t seq_len, qual_len;
    uint64_t char_len, cl, line_len, ll;
    enum read_state {OUT_READ, IN_NAME, IN_SEQ, SEQ_END, IN_QUAL} state;

    idx = (faidx_t*)calloc(1, sizeof(faidx_t));
    idx->hash = kh_init(s);
    idx->format = FAI_NONE;

    state = OUT_READ, read_done = 0, line_num = 1;
    seq_offset = qual_offset = seq_len = qual_len = char_len = cl = line_len = ll = 0;

    while ((c = bgzf_getc(bgzf)) >= 0) {
        switch (state) {
            case OUT_READ:
                switch (c) {
                    case '>':
                        if (idx->format == FAI_FASTQ) {
                            hts_log_error("Found '>' in a FASTQ file, error at line %d", line_num);
                            goto fail;
                        }

                        idx->format = FAI_FASTA;
                        state = IN_NAME;
                    break;

                    case '@':
                        if (idx->format == FAI_FASTA) {
                            hts_log_error("Found '@' in a FASTA file, error at line %d", line_num);
                            goto fail;
                        }

                        idx->format = FAI_FASTQ;
                        state = IN_NAME;
                    break;

                    case '\r':
                        // Blank line with cr-lf ending?
                        if ((c = bgzf_getc(bgzf)) == '\n') {
                            line_num++;
                        } else {
                            hts_log_error("Format error, carriage return not followed by new line at line %d", line_num);
                            goto fail;
                        }
                    break;

                    case '\n':
                        // just move onto the next line
                        line_num++;
                    break;

                    default: {
                        char s[4] = { '"', c, '"', '\0' };
                        hts_log_error("Format error, unexpected %s at line %d", isprint(c) ? s : "character", line_num);
                        goto fail;
                    break;
                    }
                }
            break;

            case IN_NAME:
                if (read_done) {
                    if (fai_insert_index(idx, name.s, seq_len, line_len, char_len, seq_offset, qual_offset) != 0)
                        goto fail;

                    read_done = 0;
                }

                name.l = 0;

                do {
                    if (!isspace(c)) {
                        kputc(c, &name);
                    } else if (name.l > 0 || c == '\n') {
                        break;
                    }
                } while ((c = bgzf_getc(bgzf)) >= 0);

                kputsn("", 0, &name);

                if (c < 0) {
                    hts_log_error("The last entry '%s' has no sequence", name.s);
                    goto fail;
                }

                // read the rest of the line if necessary
                if (c != '\n') while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');

                state = IN_SEQ; seq_len = qual_len = char_len = line_len = 0;
                seq_offset = bgzf_utell(bgzf);
                line_num++;
            break;

            case IN_SEQ:
                if (idx->format == FAI_FASTA) {
                    if (c == '\n') {
                        state = OUT_READ;
                        line_num++;
                        continue;
                    } else if (c == '>') {
                        state = IN_NAME;
                        continue;
                    }
                } else if (idx->format == FAI_FASTQ) {
                    if (c == '+') {
                        state = IN_QUAL;
                        if (c != '\n') while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');
                        qual_offset = bgzf_utell(bgzf);
                        line_num++;
                        continue;
                    } else if (c == '\n') {
                        hts_log_error("Inlined empty line is not allowed in sequence '%s' at line %d", name.s, line_num);
                        goto fail;
                    }
                }

                ll = cl = 0;

                if (idx->format == FAI_FASTA) read_done = 1;

                do {
                    ll++;
                    if (isgraph(c)) cl++;
                } while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');

                ll++; seq_len += cl;

                if (line_len == 0) {
                    line_len = ll;
                    char_len = cl;
                } else if (line_len > ll) {

                    if (idx->format == FAI_FASTA)
                        state = OUT_READ;
                    else
                        state = SEQ_END;

                } else if (line_len < ll) {
                    hts_log_error("Different line length in sequence '%s'", name.s);
                    goto fail;
                }

                line_num++;
            break;

            case SEQ_END:
                if (c == '+') {
                    state = IN_QUAL;
                    if (c != '\n') while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');
                    qual_offset = bgzf_utell(bgzf);
                    line_num++;
                    continue;
                } else {
                    hts_log_error("Format error, expecting '+', got '%c' at line %d", c, line_num);
                    goto fail;
                }
            break;

            case IN_QUAL:
                if (c == '\n') {
                    if (!read_done) {
                        hts_log_error("Inlined empty line is not allowed in quality of sequence '%s'", name.s);
                        goto fail;
                    }

                    state = OUT_READ;
                    line_num++;
                    continue;
                } else if (c == '@' && read_done) {
                    state = IN_NAME;
                    continue;
                }

                ll = cl = 0;

                do {
                    ll++;
                    if (isgraph(c)) cl++;
                } while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');

                ll++; qual_len += cl;

                if (line_len < ll) {
                    hts_log_error("Quality line length too long in '%s' at line %d", name.s, line_num);
                    goto fail;
                } else if (qual_len == seq_len) {
                    read_done = 1;
                } else if (qual_len > seq_len) {
                    hts_log_error("Quality length longer than sequence in '%s' at line %d", name.s, line_num);
                    goto fail;
                } else if (line_len > ll) {
                    hts_log_error("Quality line length too short in '%s' at line %d", name.s, line_num);
                    goto fail;
                }

                line_num++;
            break;
        }
    }

    if (read_done) {
        if (fai_insert_index(idx, name.s, seq_len, line_len, char_len, seq_offset, qual_offset) != 0)
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

        if (fai->format == FAI_FASTA) {
            snprintf(buf, sizeof(buf),
                 "\t%"PRIu64"\t%"PRIu64"\t%"PRIu32"\t%"PRIu32"\n",
                 x.len, x.seq_offset, x.line_blen, x.line_len);
        } else {
            snprintf(buf, sizeof(buf),
                 "\t%"PRIu64"\t%"PRIu64"\t%"PRIu32"\t%"PRIu32"\t%"PRIu64"\n",
                 x.len, x.seq_offset, x.line_blen, x.line_len, x.qual_offset);
        }

        if (hputs(fai->name[i], fp) != 0) return -1;
        if (hputs(buf, fp) != 0) return -1;
    }
    return 0;
}


static faidx_t *fai_read(hFILE *fp, const char *fname, int format)
{
    faidx_t *fai;
    char *buf = NULL, *p;
    ssize_t l, lnum = 1;

    fai = (faidx_t*)calloc(1, sizeof(faidx_t));
    if (!fai) return NULL;

    fai->hash = kh_init(s);
    if (!fai->hash) goto fail;

    buf = (char*)calloc(0x10000, 1);
    if (!buf) goto fail;

    while ((l = hgetln(buf, 0x10000, fp)) > 0) {
        uint32_t line_len, line_blen, n;
        uint64_t len;
        uint64_t seq_offset;
        uint64_t qual_offset = 0;

        for (p = buf; *p && !isspace_c(*p); ++p);

        if (p - buf < l) {
            *p = 0; ++p;
        }

        if (format == FAI_FASTA) {
            n = sscanf(p, "%"SCNu64"%"SCNu64"%"SCNu32"%"SCNu32"", &len, &seq_offset, &line_blen, &line_len);

            if (n != 4) {
                hts_log_error("Could not understand FASTA index %s line %zd", fname, lnum);
                goto fail;
            }
        } else {
            n = sscanf(p, "%"SCNu64"%"SCNu64"%"SCNu32"%"SCNu32"%"SCNu64"", &len, &seq_offset, &line_blen, &line_len, &qual_offset);

            if (n != 5) {
                if (n == 4) {
                    hts_log_error("Possibly this is a FASTA index, try using faidx.  Problem in %s line %zd", fname, lnum);
                } else {
                    hts_log_error("Could not understand FASTQ index %s line %zd", fname, lnum);
                }

                goto fail;
            }
        }

        if (fai_insert_index(fai, buf, len, line_len, line_blen, seq_offset, qual_offset) != 0) {
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


static int fai_build3_core(const char *fn, const char *fnfai, const char *fngzi)
{
    kstring_t fai_kstr = { 0, 0, NULL };
    kstring_t gzi_kstr = { 0, 0, NULL };
    BGZF *bgzf = NULL;
    hFILE *fp = NULL;
    faidx_t *fai = NULL;
    int save_errno, res;
    char *file_type;

    bgzf = bgzf_open(fn, "r");

    if ( !bgzf ) {
        hts_log_error("Failed to open the file %s", fn);
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

    if (fai->format == FAI_FASTA) {
        file_type   = "FASTA";
    } else {
        file_type   = "FASTQ";
    }

    if (!fnfai) {
        if (ksprintf(&fai_kstr, "%s.fai", fn) < 0) goto fail;
        fnfai = fai_kstr.s;
    }

    if (!fngzi) {
        if (ksprintf(&gzi_kstr, "%s.gzi", fn) < 0) goto fail;
        fngzi = gzi_kstr.s;
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
        hts_log_error("Failed to open %s index %s : %s", file_type, fnfai, strerror(errno));
        goto fail;
    }

    if (fai_save(fai, fp) != 0) {
        hts_log_error("Failed to write %s index %s : %s", file_type, fnfai, strerror(errno));
        goto fail;
    }

    if (hclose(fp) != 0) {
        hts_log_error("Failed on closing %s index %s : %s", file_type, fnfai, strerror(errno));
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


int fai_build3(const char *fn, const char *fnfai, const char *fngzi) {
    return fai_build3_core(fn, fnfai, fngzi);
}


int fai_build(const char *fn) {
    return fai_build3(fn, NULL, NULL);
}


static faidx_t *fai_load3_core(const char *fn, const char *fnfai, const char *fngzi,
                   int flags, int format)
{
    kstring_t fai_kstr = { 0, 0, NULL };
    kstring_t gzi_kstr = { 0, 0, NULL };
    hFILE *fp = NULL;
    faidx_t *fai = NULL;
    int res, gzi_index_needed = 0;
    char *file_type;

    if (format == FAI_FASTA) {
        file_type   = "FASTA";
    } else {
        file_type   = "FASTQ";
    }

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

    if (fp) {
        // index file present, check if a compressed index is needed
        hFILE *gz = NULL;
        BGZF *bgzf = bgzf_open(fn, "rb");

        if (bgzf == 0) {
            hts_log_error("Failed to open %s file %s", file_type, fn);
            goto fail;
        }

        if (bgzf_compression(bgzf) == 2) { // BGZF compression
            if ((gz = hopen(fngzi, "rb")) == 0) {

                if (!(flags & FAI_CREATE) || errno != ENOENT) {
                    hts_log_error("Failed to open %s index %s: %s", file_type, fngzi, strerror(errno));
                    bgzf_close(bgzf);
                    goto fail;
                }

                gzi_index_needed = 1;
                res = hclose(fp); // closed as going to be re-indexed

                if (res < 0) {
                    hts_log_error("Failed on closing %s index %s : %s", file_type, fnfai, strerror(errno));
                    goto fail;
                }
            } else {
                res = hclose(gz);

                if (res < 0) {
                    hts_log_error("Failed on closing %s index %s : %s", file_type, fngzi, strerror(errno));
                    goto fail;
                }
            }
        }

        bgzf_close(bgzf);
    }

    if (fp == 0 || gzi_index_needed) {
        if (!(flags & FAI_CREATE) || errno != ENOENT) {
            hts_log_error("Failed to open %s index %s: %s", file_type, fnfai, strerror(errno));
            goto fail;
        }

        hts_log_info("Build %s index", file_type);

        if (fai_build3_core(fn, fnfai, fngzi) < 0) {
            goto fail;
        }

        fp = hopen(fnfai, "rb");
        if (fp == 0) {
            hts_log_error("Failed to open %s index %s: %s", file_type, fnfai, strerror(errno));
            goto fail;
        }
    }

    fai = fai_read(fp, fnfai, format);
    if (fai == NULL) {
        hts_log_error("Failed to read %s index %s", file_type, fnfai);
        goto fail;
    }

    res = hclose(fp);
    fp = NULL;
    if (res < 0) {
        hts_log_error("Failed on closing %s index %s : %s", file_type, fnfai, strerror(errno));
        goto fail;
    }

    fai->bgzf = bgzf_open(fn, "rb");
    if (fai->bgzf == 0) {
        hts_log_error("Failed to open %s file %s", file_type, fn);
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


faidx_t *fai_load3(const char *fn, const char *fnfai, const char *fngzi,
                   int flags) {
    return fai_load3_core(fn, fnfai, fngzi, flags, FAI_FASTA);
}


faidx_t *fai_load(const char *fn)
{
    return fai_load3(fn, NULL, NULL, FAI_CREATE);
}


faidx_t *fai_load3_format(const char *fn, const char *fnfai, const char *fngzi,
                   int flags, enum fai_format_options format) {
    return fai_load3_core(fn, fnfai, fngzi, flags, format);
}


faidx_t *fai_load_format(const char *fn, enum fai_format_options format) {
    return fai_load3_format(fn, NULL, NULL, FAI_CREATE, format);
}


static char *fai_retrieve(const faidx_t *fai, const faidx1_t *val,
                          uint64_t offset, long beg, long end, int *len) {
    char *s;
    size_t l;
    int c = 0;
    int ret = bgzf_useek(fai->bgzf,
                         offset
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


static int fai_get_val(const faidx_t *fai, const char *str, int *len, faidx1_t *val, long *fbeg, long *fend) {
    char *s, *ep;
    size_t i, l, k, name_end;
    khiter_t iter;
    khash_t(s) *h;
    long beg, end;

    beg = end = -1;
    h = fai->hash;
    name_end = l = strlen(str);
    s = (char*)malloc(l+1);
    if (!s) {
        *len = -1;
        return 1;
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
        hts_log_warning("Reference %s not found in file, returning empty sequence", str);
        free(s);
        *len = -2;
        return 1;
    }
    *val = kh_value(h, iter);
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
        end = i < k? strtol(s + i, &ep, 10) : val->len;
        if (beg > 0) --beg;
        // Check for out of range numbers.  Only going to be a problem on
        // 32-bit platforms with >2Gb sequence length.
        if (errno == ERANGE && (uint64_t) val->len > LONG_MAX) {
            hts_log_error("Positions in range %s are too large for this platform", s);
            free(s);
            *len = -3;
            return 1;
        }
        errno = save_errno;
    } else beg = 0, end = val->len;
    if (beg >= val->len) beg = val->len;
    if (end >= val->len) end = val->len;
    if (beg > end) beg = end;
    free(s);

    *fbeg = beg;
    *fend = end;

    return 0;
}


char *fai_fetch(const faidx_t *fai, const char *str, int *len)
{
    faidx1_t val;
    long beg, end;

    if (fai_get_val(fai, str, len, &val, &beg, &end)) {
        return NULL;
    }

    // now retrieve the sequence
    return fai_retrieve(fai, &val, val.seq_offset, beg, end, len);
}


char *fai_fetchqual(const faidx_t *fai, const char *str, int *len) {
    faidx1_t val;
    long beg, end;

    if (fai_get_val(fai, str, len, &val, &beg, &end)) {
        return NULL;
    }

    // now retrieve the sequence
    return fai_retrieve(fai, &val, val.qual_offset, beg, end, len);
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


static int faidx_adjust_position(const faidx_t *fai, faidx1_t *val, const char *c_name, int *p_beg_i, int *p_end_i, int *len) {
    khiter_t iter;

    // Adjust position
    iter = kh_get(s, fai->hash, c_name);

    if (iter == kh_end(fai->hash)) {
        *len = -2;
        hts_log_error("The sequence \"%s\" not found", c_name);
        return 1;
    }

    *val = kh_value(fai->hash, iter);

    if(*p_end_i < *p_beg_i)
        *p_beg_i = *p_end_i;

    if(*p_beg_i < 0)
        *p_beg_i = 0;
    else if(val->len <= *p_beg_i)
        *p_beg_i = val->len - 1;

    if(*p_end_i < 0)
        *p_end_i = 0;
    else if(val->len <= *p_end_i)
        *p_end_i = val->len - 1;

    return 0;
}


char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
    faidx1_t val;

    // Adjust position
    if (faidx_adjust_position(fai, &val, c_name, &p_beg_i, &p_end_i, len)) {
        return NULL;
    }

    // Now retrieve the sequence
    return fai_retrieve(fai, &val, val.seq_offset, p_beg_i, (long) p_end_i + 1, len);
}


char *faidx_fetch_qual(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
    faidx1_t val;

    // Adjust position
    if (faidx_adjust_position(fai, &val, c_name, &p_beg_i, &p_end_i, len)) {
        return NULL;
    }

    // Now retrieve the sequence
    return fai_retrieve(fai, &val, val.qual_offset, p_beg_i, (long) p_end_i + 1, len);
}


int faidx_has_seq(const faidx_t *fai, const char *seq)
{
    khiter_t iter = kh_get(s, fai->hash, seq);
    if (iter == kh_end(fai->hash)) return 0;
    return 1;
}

