/* test/test_bgzf.c -- bgzf unit tests

   Copyright (C) 2017 Genome Research Ltd

   Author: Robert Davies <rmd@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE.
 */

#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "hfile_internal.h"

const char *bgzf_suffix = ".gz";
const char *idx_suffix  = ".gzi";
const char *tmp_suffix  = ".tmp";

#define BUFSZ 32768

typedef struct {
    char *src_plain;
    char *src_bgzf;
    char *src_idx;
    char *tmp_bgzf;
    char *tmp_idx;
    FILE *f_plain;
    FILE *f_bgzf;
    FILE *f_idx;
    const unsigned char *text;
    size_t ltext;
} Files;

typedef enum {
    USE_BGZF_OPEN,
    USE_BGZF_DOPEN,
    USE_BGZF_HOPEN
} Open_method;

static FILE * try_fopen(const char *name, const char *mode) {
    FILE *f = fopen(name, mode);
    if (!f) {
        fprintf(stderr, "Couldn't open %s : %s\n", name, strerror(errno));
        return NULL;
    }
    return f;
}

static int try_fclose(FILE **file, const char *name, const char *func) {
    FILE *to_close = *file;
    *file = NULL;
    if (fclose(to_close) != 0) {
        fprintf(stderr, "%s : Error on closing %s : %s\n",
                func, name, strerror(errno));
        return -1;
    }

    return 0;
}

static ssize_t try_fread(FILE *in, void *buf, size_t len,
                         const char *func, const char *fname) {
    size_t got = fread(buf, 1, len, in);
    if (got == 0 && ferror(in)) {
        fprintf(stderr, "%s : Error reading from %s : %s\n",
                func, fname, strerror(errno));
        return -1;
    }
    return got;
}

static int try_fseek_start(FILE *f, const char *name, const char *func) {
    if (0 != fseek(f, 0, SEEK_SET)) {
        fprintf(stderr, "%s : Couldn't seek on %s : %s\n",
                func, name, strerror(errno));
        return -1;
    }
    return 0;
}

static BGZF * try_bgzf_open(const char *name, const char *mode,
                            const char *func) {
    BGZF * bgz = bgzf_open(name, mode);
    if (!bgz) {
        fprintf(stderr, "%s : Couldn't bgzf_open %s with mode %s : %s\n",
                func, name, mode, strerror(errno));
        return NULL;
    }
    return bgz;
}

static BGZF * try_bgzf_dopen(const char *name, const char *mode,
                             const char *func) {
    BGZF *bgz = NULL;
    int fd = open(name, hfile_oflags(mode), 0666);
    if (fd < 0) {
        fprintf(stderr, "%s : Failed to open %s with mode %s : %s\n",
                func, name, mode, strerror(errno));
        return NULL;
    }

    bgz = bgzf_dopen(fd, mode);
    if (!bgz) {
        fprintf(stderr, "%s : bgzf_dopen failed on %s mode %s : %s\n",
                func, name, mode, strerror(errno));
        close(fd);
        return NULL;
    }

    return bgz;
}

static BGZF * try_bgzf_hopen(const char *name, const char *mode,
                          const char *func) {
    hFILE *hfp = hopen(name, mode);
    BGZF *bgz = NULL;

    if (!hfp) {
        fprintf(stderr, "%s : hopen failed on %s mode %s : %s\n",
                func, name, mode, strerror(errno));
        return NULL;
    }

    bgz = bgzf_hopen(hfp, mode);
    if (!bgz) {
        fprintf(stderr, "%s : bgzf_hopen failed on %s mode %s : %s\n",
                func, name, mode, strerror(errno));
        hclose_abruptly(hfp);
        return NULL;
    }

    return bgz;
}

static int try_bgzf_close(BGZF **bgz, const char *name, const char *func) {
    BGZF *to_close = *bgz;
    *bgz = NULL;
    if (bgzf_close(to_close) != 0) {
        fprintf(stderr, "%s : bgzf_close failed on %s : %s\n",
                func, name, strerror(errno));
        return -1;
    }
    return 0;
}

static ssize_t try_bgzf_read(BGZF *fp, void *data, size_t length,
                             const char *name, const char *func) {
    ssize_t got = bgzf_read(fp, data, length);
    if (got < 0) {
        fprintf(stderr, "%s : Error from bgzf_read %s : %s\n",
                func, name, strerror(errno));
    }
    return got;
}

static ssize_t try_bgzf_write(BGZF *fp, const void *data, size_t length,
                              const char *name, const char *func) {
    ssize_t put = bgzf_write(fp, data, length);
    if (put < (ssize_t) length) {
        fprintf(stderr, "%s : %s %s : %s\n",
                func, put < 0 ? "Error writing to" : "Short write on",
                name, strerror(errno));
        return -1;
    }

    return put;
}

static int try_bgzf_compression(BGZF *fp, int expect,
                                const char *name, const char *func) {
    int res = bgzf_compression(fp);
    if (res != expect) {
        fprintf(stderr,
                "%s : Unexpected result %d from bgzf_compression on %s; "
                "expected %d\n",
                func, res, name, expect);
        return -1;
    }
    return 0;
}

static int try_bgzf_mt(BGZF *bgz, int nthreads, const char *func) {
    if (bgzf_mt(bgz, nthreads, 64) != 0) {
        fprintf(stderr, "%s : Error from bgzf_mt : %s\n",
                func, strerror(errno));
        return -1;
    }
    return 0;
}

static int try_bgzf_index_build_init(BGZF *bgz,
                                     const char *name, const char *func) {
    if (bgzf_index_build_init(bgz) != 0) {
        fprintf(stderr, "%s : Error from bgzf_index_build_init on %s : %s\n",
                func, name, strerror(errno));
        return -1;
    }
    return 0;
}

static int try_bgzf_index_load(BGZF *fp, const char *bname, const char *suffix,
                               const char *func) {
    if (bgzf_index_load(fp, bname, suffix) != 0) {
        fprintf(stderr, "%s : Couldn't bgzf_index_load %s%s : %s\n",
                func, bname, suffix ? suffix : "", strerror(errno));
        return -1;
    }
    return 0;
}

static int try_bgzf_index_dump(BGZF *fp, const char *bname, const char *suffix,
                               const char *func) {
    if (bgzf_index_dump(fp, bname, suffix) != 0) {
        fprintf(stderr, "%s : Couldn't bgzf_index_dump %s%s : %s\n",
                func, bname, suffix ? suffix : "", strerror(errno));
        return -1;
    }
    return 0;
}

static int try_bgzf_useek(BGZF *fp, long uoffset, int where,
                          const char *name, const char *func) {
    if (bgzf_useek(fp, uoffset, where) < 0) {
        fprintf(stderr, "%s : Error from bgzf_useek(%s, %ld, %d) : %s\n",
                func, name, uoffset, where, strerror(errno));
        return -1;
    }
    return 0;
}

static int try_bgzf_getc(BGZF *fp, size_t pos, int expected,
                         const char *name, const char *func) {
    int c = bgzf_getc(fp);
    if (c != expected) {
        fprintf(stderr,
                "%s : Unexpected value (%d) from bgzf_getc on %s pos %zu; "
                "expected %d\n",
                func, c, name, pos, expected);
        return -1;
    }
    return c;
}

static int compare_buffers(const unsigned char *b1, const unsigned char *b2,
                           size_t l1, size_t l2,
                           const char *name1, const char *name2,
                           const char *func) {
    if (l1 != l2) {
        fprintf(stderr, "%s : EOF on %s\n", func, l1 < l2 ? name1 : name2);
        return -1;
    }
    if (memcmp(b1, b2, l1) != 0) {
        fprintf(stderr, "%s : difference between %s and %s\n",
                func, name1, name2);
        return -1;
    }

    return 0;
}

static void cleanup(Files *f, int retval) {
    /* Remove temp files if successful.  If not, leave them for inspection */
    if (retval == EXIT_SUCCESS) {
        unlink(f->tmp_bgzf);
        unlink(f->tmp_idx);
    }
    if (f->f_plain) fclose(f->f_plain);
    if (f->f_bgzf)  fclose(f->f_bgzf);
    if (f->f_idx)   fclose(f->f_idx);
    free(f->src_plain);
    free((unsigned char *) f->text);
}

static int setup(const char *src, Files *f) {
    size_t len = (strlen(src) + strlen(bgzf_suffix) + strlen(idx_suffix)
                  + strlen(tmp_suffix) + 8);
    char *mem, *text;
    const unsigned int max = 50000;
    unsigned int i;
    size_t text_sz = max * 8 + 1;

    mem = calloc(5, len);
    if (mem == NULL) {
        perror(__func__);
        return -1;
    }

    snprintf(mem,           len, "%s",     src);
    snprintf(mem + len * 1, len, "%s%s",   src, bgzf_suffix);
    snprintf(mem + len * 2, len, "%s%s%s", src, bgzf_suffix, idx_suffix);
    snprintf(mem + len * 3, len, "%s%s%s", src, tmp_suffix, bgzf_suffix);
    snprintf(mem + len * 4, len, "%s%s%s%s",
             src, tmp_suffix, bgzf_suffix, idx_suffix);

    f->src_plain  = mem;
    f->src_bgzf   = mem + len * 1;
    f->src_idx    = mem + len * 2;
    f->tmp_bgzf  = mem + len * 3;
    f->tmp_idx   = mem + len * 4;

    text = malloc(text_sz);
    if (!text) {
        perror(__func__);
        goto fail;
    }
    for (i = 0; i < max; i++) snprintf(text + i*8, text_sz - i*8, "%07d\n", i);
    f->text = (unsigned char *) text;
    f->ltext = text_sz - 1;

    if ((f->f_plain = try_fopen(f->src_plain, "rb")) == NULL) goto fail;
    if ((f->f_bgzf  = try_fopen(f->src_bgzf,  "rb")) == NULL) goto fail;
    if ((f->f_idx   = try_fopen(f->src_idx,   "rb")) == NULL) goto fail;

    return 0;

 fail:
    return -1;
}

static int test_read(Files *f) {
    BGZF* bgz;
    ssize_t bg_got, f_got;
    unsigned char bg_buf[BUFSZ], f_buf[BUFSZ];

    bgz = try_bgzf_open(f->src_bgzf, "r", __func__);
    if (!bgz) return -1;

    do {
        bg_got = try_bgzf_read(bgz, bg_buf, BUFSZ, f->src_bgzf, __func__);
        if (bg_got < 0) goto fail;

        f_got = try_fread(f->f_plain, f_buf, BUFSZ, __func__, f->src_plain);
        if (f_got < 0) goto fail;

        if (compare_buffers(f_buf, bg_buf, f_got, bg_got,
                            f->src_plain, f->src_bgzf, __func__) != 0) {
            goto fail;
        }
    } while (bg_got > 0 && f_got > 0);

    if (try_bgzf_close(&bgz, f->src_bgzf, __func__) != 0) return -1;
    if (try_fseek_start(f->f_plain, f->src_plain, __func__) != 0) return -1;

    return 0;

 fail:
    if (bgz) bgzf_close(bgz);
    return -1;
}

static int test_write_read(Files *f, const char *mode, Open_method method,
                           int nthreads, int expected_compression) {
    BGZF* bgz = NULL;
    ssize_t bg_put, bg_got;
    size_t pos = 0;
    unsigned char bg_buf[BUFSZ];

    switch (method) {
    case USE_BGZF_DOPEN:
        bgz = try_bgzf_dopen(f->tmp_bgzf, mode, __func__);
        break;
    case USE_BGZF_HOPEN:
        bgz = try_bgzf_hopen(f->tmp_bgzf, mode, __func__);
        break;
    default:
        bgz = try_bgzf_open(f->tmp_bgzf, mode, __func__);
        break;
    }
    if (!bgz) goto fail;

    if (nthreads > 0 && try_bgzf_mt(bgz, nthreads, __func__) != 0) goto fail;

    bg_put = try_bgzf_write(bgz, f->text, f->ltext, f->tmp_bgzf, __func__);
    if (bg_put < 0) goto fail;

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;

    switch (method) {
    case USE_BGZF_DOPEN:
        bgz = try_bgzf_dopen(f->tmp_bgzf, "r", __func__);
        break;
    case USE_BGZF_HOPEN:
        bgz = try_bgzf_hopen(f->tmp_bgzf, "r", __func__);
        break;
    default:
        bgz = try_bgzf_open(f->tmp_bgzf, "r", __func__);
        break;
    }
    if (!bgz) goto fail;

    if (nthreads > 0 && try_bgzf_mt(bgz, nthreads, __func__) != 0) goto fail;

    if (try_bgzf_compression(bgz, expected_compression,
                             f->tmp_bgzf, __func__) != 0) {
        goto fail;
    }

    do {
        bg_got = try_bgzf_read(bgz, bg_buf, BUFSZ, f->tmp_bgzf, __func__);
        if (bg_got < 0) goto fail;

        if (pos < f->ltext &&
            memcmp(f->text + pos, bg_buf,
                   pos + bg_got < f->ltext ? bg_got : f->ltext - pos) != 0) {
            fprintf(stderr, "%s : Got wrong data from %s, pos %zu\n",
                    __func__, f->tmp_bgzf, pos);
            goto fail;
        }
        pos += bg_got;
    } while (bg_got > 0);

    if (pos != bg_put) {
        fprintf(stderr, "%s : bgzf_read got %zd bytes; expected %zd\n",
                __func__, pos, bg_put);
        goto fail;
    }

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;

    return 0;

 fail:
    if (bgz) bgzf_close(bgz);
    return -1;
}

static int test_embed_eof(Files *f, const char *mode, int nthreads) {
    BGZF* bgz = NULL;
    ssize_t bg_put, bg_got;
    size_t pos = 0, half = BUFSZ < f->ltext ? BUFSZ : f->ltext / 2;
    char append_mode[16];
    unsigned char bg_buf[BUFSZ];

    for (pos = 0; pos < sizeof(append_mode) - 1 && mode[pos] != 0; pos++) {
        append_mode[pos] = mode[pos] == 'w' ? 'a' : mode[pos];
    }
    append_mode[pos] ='\0';

    // Write first half
    bgz = try_bgzf_open(f->tmp_bgzf, mode, __func__);
    if (!bgz) goto fail;

    if (nthreads > 0 && try_bgzf_mt(bgz, nthreads, __func__) != 0) goto fail;

    bg_put = try_bgzf_write(bgz, f->text, half, f->tmp_bgzf, __func__);
    if (bg_put < 0) goto fail;

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;


    // Write second half.  Append mode, so an EOF block should be in the
    // middle.
    bgz = try_bgzf_open(f->tmp_bgzf, append_mode, __func__);
    if (!bgz) goto fail;

    if (nthreads > 0 && try_bgzf_mt(bgz, nthreads, __func__) != 0) goto fail;

    bg_put = try_bgzf_write(bgz, f->text + half, f->ltext - half, f->tmp_bgzf,
                            __func__);
    if (bg_put < 0) goto fail;

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;

    // Try reading
    pos = 0;
    bgz = try_bgzf_open(f->tmp_bgzf, "r", __func__);
    if (!bgz) goto fail;

    if (nthreads > 0 && try_bgzf_mt(bgz, nthreads, __func__) != 0) goto fail;

    do {
        bg_got = try_bgzf_read(bgz, bg_buf, BUFSZ, f->tmp_bgzf, __func__);
        if (bg_got < 0) goto fail;

        if (pos < f->ltext &&
            memcmp(f->text + pos, bg_buf,
                   pos + bg_got < f->ltext ? bg_got : f->ltext - pos) != 0) {
            fprintf(stderr, "%s : Got wrong data from %s, pos %zu\n",
                    __func__, f->tmp_bgzf, pos);
            goto fail;
        }
        pos += bg_got;
    } while (bg_got > 0);

    if (pos != f->ltext) {
        fprintf(stderr, "%s : bgzf_read got %zd bytes; expected %zd\n",
                __func__, pos, f->ltext);
        goto fail;
    }

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;

    return 0;

 fail:
    if (bgz) bgzf_close(bgz);
    return -1;
}

static int test_index_load_dump(Files *f) {
    BGZF* bgz = NULL;
    FILE *fdest = NULL;
    unsigned char buf_src[BUFSZ], buf_dest[BUFSZ];
    ssize_t got_src, got_dest;

    bgz = try_bgzf_open(f->src_bgzf, "r", __func__);
    if (!bgz) return -1;

    if (try_bgzf_index_load(bgz, f->src_bgzf, idx_suffix, __func__) != 0) {
        goto fail;
    }

    if (try_bgzf_index_dump(bgz, f->tmp_bgzf, idx_suffix, __func__) != 0) {
        goto fail;
    }

    fdest = try_fopen(f->tmp_idx, "r");
    do {
        got_src  = try_fread(f->f_idx, buf_src,  BUFSZ, __func__, f->src_idx);
        if (got_src < 0) goto fail;
        got_dest = try_fread(fdest,    buf_dest, BUFSZ, __func__, f->tmp_idx);
        if (got_dest < 0) goto fail;
        if (compare_buffers(buf_src, buf_dest, got_src, got_dest,
                            f->src_idx, f->tmp_idx, __func__) != 0) goto fail;
    } while (got_src > 0 && got_dest > 0);
    if (try_fclose(&fdest, f->tmp_idx, __func__) != 0) goto fail;

    if (try_bgzf_close(&bgz, f->src_bgzf, __func__) != 0) goto fail;

    return 0;

 fail:
    if (fdest) fclose(fdest);
    if (bgz) bgzf_close(bgz);
    return -1;
}

static int test_check_EOF(char *name, int expected) {
    BGZF *bgz = try_bgzf_open(name, "r", __func__);
    int eof;
    if (!bgz) return -1;
    eof = bgzf_check_EOF(bgz);
    if (eof != expected) {
        fprintf(stderr, "%s : Unexpected result %d from bgzf_check_EOF on %s; "
                "expected %d\n",
                __func__, eof, name, expected);
        bgzf_close(bgz);
        return -1;
    }

    return try_bgzf_close(&bgz, name, __func__);
}

static int test_index_seek_getc(Files *f, const char *mode,
                                int cache_size, int nthreads) {
    BGZF* bgz = NULL;
    ssize_t bg_put;
    size_t i, j, iskip = f->ltext / 10;
    int is_uncompressed = strchr(mode, 'u') != NULL;

    bgz = try_bgzf_open(f->tmp_bgzf, mode, __func__);
    if (!bgz) goto fail;

    if (try_bgzf_index_build_init(bgz, f->tmp_bgzf, __func__) != 0) goto fail;

    if (nthreads > 0 && try_bgzf_mt(bgz, nthreads, __func__) != 0) goto fail;

    bg_put = try_bgzf_write(bgz, f->text, f->ltext, f->tmp_bgzf, __func__);
    if (bg_put < 0) goto fail;

    if (!is_uncompressed) {
        if (try_bgzf_index_dump(bgz, f->tmp_idx, NULL, __func__) != 0) {
            goto fail;
        }
    }

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;

    bgz = try_bgzf_open(f->tmp_bgzf, "r", __func__);
    if (!bgz) goto fail;

    if (nthreads > 0 && try_bgzf_mt(bgz, nthreads, __func__) != 0) goto fail;

    if (!is_uncompressed) {
        if (try_bgzf_index_load(bgz, f->tmp_bgzf, idx_suffix, __func__) != 0) {
            goto fail;
        }
    }

    for (i = 0; i < f->ltext; i += iskip) {
        if (try_bgzf_useek(bgz, i, SEEK_SET, f->tmp_bgzf, __func__) != 0) {
            goto fail;
        }

        for (j = 0; j < 16 && i + j < f->ltext; j++) {
            if (try_bgzf_getc(bgz, i + j, f->text[i + j],
                              f->tmp_bgzf, __func__) < 0) goto fail;
        }
    }

    if (try_bgzf_useek(bgz, 0, SEEK_SET, f->tmp_bgzf, __func__) != 0) {
        goto fail;
    }
    for (j = 0; j < 70000 && j < f->ltext; j++) { // Should force a block load
        if (try_bgzf_getc(bgz, j, f->text[j],
                          f->tmp_bgzf, __func__) < 0) goto fail;
    }

    if (cache_size > 0) {
        size_t mid = f->ltext / 2;
        bgzf_set_cache_size(bgz, cache_size);

        for (i = 0; i < 10; i++) {
            if (try_bgzf_useek(bgz, 0, SEEK_SET, f->tmp_bgzf, __func__) != 0) {
                goto fail;
            }
            for (j = 0; j < 64 && j < f->ltext; j++) {
                if (try_bgzf_getc(bgz, j, f->text[j],
                                  f->tmp_bgzf, __func__) < 0) goto fail;
            }

            if (try_bgzf_useek(bgz, mid, SEEK_SET,
                               f->tmp_bgzf, __func__) != 0) {
                goto fail;
            }
            for (j = 0; j < 64 && j + mid < f->ltext; j++) {
                if (try_bgzf_getc(bgz, j + mid, f->text[j + mid],
                                  f->tmp_bgzf, __func__) < 0) goto fail;
            }
        }
    }

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;

    return 0;

 fail:
    if (bgz) bgzf_close(bgz);
    return -1;
}

static int test_bgzf_getline(Files *f, const char *mode, int nthreads) {
    BGZF* bgz = NULL;
    ssize_t bg_put;
    size_t pos;
    kstring_t str = { 0, 0, NULL };
    const char *text = (const char *) f->text;

    bgz = try_bgzf_open(f->tmp_bgzf, mode, __func__);
    if (!bgz) goto fail;

    if (nthreads > 0 && try_bgzf_mt(bgz, nthreads, __func__) != 0) goto fail;

    bg_put = try_bgzf_write(bgz, f->text, f->ltext, f->tmp_bgzf, __func__);
    if (bg_put < 0) goto fail;

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;

    bgz = try_bgzf_open(f->tmp_bgzf, "r", __func__);
    if (!bgz) goto fail;

    for (pos = 0; pos < f->ltext; ) {
        const char *end = strchr(text + pos, '\n');
        size_t l = end ? end - (text + pos) : f->ltext - pos;
        int res;

        if ((res = bgzf_getline(bgz, '\n', &str)) < 0) {
            fprintf(stderr, "%s : %s from bgzf_getline on %s : %s\n",
                    __func__, res < -1 ? "Error" : "Unexpected EOF",
                    f->tmp_bgzf, res < -1 ? strerror(errno) : "EOF");
            goto fail;
        }

        if (str.l != l || memcmp(text + pos, str.s, l) != 0) {
            fprintf(stderr,
                    "%s : Unexpected data from bgzf_getline on %s\n"
                    "Expected : %.*s\n"
                    "Got      : %.*s\n",
                    __func__, f->tmp_bgzf, (int) l, (char *) f->text + pos,
                    (int) str.l, str.s);
        }

        pos += l + 1;
    }

    if (try_bgzf_close(&bgz, f->tmp_bgzf, __func__) != 0) goto fail;
    return 0;

 fail:
    if (bgz) bgzf_close(bgz);
    return -1;
}

int main(int argc, char **argv) {
    Files f = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0 };
    int retval = EXIT_FAILURE;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <source file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (setup(argv[1], &f) != 0) goto out;

    // Try reading an existing file
    if (test_check_EOF(f.src_bgzf, 1) != 0) goto out;
    if (test_read(&f) != 0) goto out;

    // Try writing some data and reading it back
    if (test_write_read(&f, "wu", USE_BGZF_OPEN, 0, 0) != 0) goto out;
    if (test_check_EOF(f.tmp_bgzf, 0) != 0) goto out;
    if (test_write_read(&f, "w",  USE_BGZF_OPEN, 0, 2) != 0)  goto out;
    if (test_check_EOF(f.tmp_bgzf, 1) != 0) goto out;
    if (test_write_read(&f, "w0", USE_BGZF_OPEN, 0, 2) != 0) goto out;
    if (test_check_EOF(f.tmp_bgzf, 1) != 0) goto out;
    if (test_write_read(&f, "w1", USE_BGZF_DOPEN, 0, 2) != 0) goto out;
    if (test_check_EOF(f.tmp_bgzf, 1) != 0) goto out;
    if (test_write_read(&f, "w9", USE_BGZF_HOPEN, 0, 2) != 0) goto out;
    if (test_check_EOF(f.tmp_bgzf, 1) != 0) goto out;
    if (test_write_read(&f, "wg", USE_BGZF_OPEN, 0, 1) != 0) goto out;
    if (test_check_EOF(f.tmp_bgzf, 0) != 0) goto out;

    // Try writing and reading with threads
    if (test_write_read(&f, "w", USE_BGZF_OPEN, 1, 2) != 0) goto out;
    if (test_check_EOF(f.tmp_bgzf, 1) != 0) goto out;
    if (test_write_read(&f, "w", USE_BGZF_OPEN, 2, 2) != 0) goto out;
    if (test_check_EOF(f.tmp_bgzf, 1) != 0) goto out;

    // Embedded EOF block
    if (test_embed_eof(&f, "w", 0) != 0) goto out;
    if (test_embed_eof(&f, "w", 1) != 0) goto out;
    if (test_embed_eof(&f, "w", 2) != 0) goto out;

    // Index load and dump
    if (test_index_load_dump(&f) != 0) goto out;

    // Index building on the fly and bgzf_useek
    if (test_index_seek_getc(&f, "w", 1000000, 0) != 0) goto out;

    // Index building on the fly and bgzf_useek, with threads
    // ** Not implemented yet **
    // if (test_index_seek_getc(&f, "w", 1000000, 1) != 0) goto out;
    // if (test_index_seek_getc(&f, "w", 1000000, 2) != 0) goto out;

    // bgzf_useek on an uncompressed file
    if (test_index_seek_getc(&f, "wu", 0, 0) != 0) goto out;

    // getline
    if (test_bgzf_getline(&f, "w", 0) != 0) goto out;
    if (test_bgzf_getline(&f, "w", 1) != 0) goto out;
    if (test_bgzf_getline(&f, "w", 2) != 0) goto out;

    retval = EXIT_SUCCESS;

 out:
    cleanup(&f, retval);
    return retval;
}
