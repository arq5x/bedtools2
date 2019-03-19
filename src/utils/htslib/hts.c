/*  hts.c -- format-neutral I/O, indexing, and iterator API functions.

    Copyright (C) 2008, 2009, 2012-2017 Genome Research Ltd.
    Copyright (C) 2012, 2013 Broad Institute.

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

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <limits.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include <assert.h>

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "cram/cram.h"
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"
#include "version.h"
#include "hts_internal.h"
#include "hfile_internal.h"
#include "htslib/hts_os.h" // drand48

#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/ksort.h"

KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

int hts_verbose = HTS_LOG_WARNING;

const char *hts_version()
{
    return HTS_VERSION;
}

const unsigned char seq_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

const int seq_nt16_int[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

/**********************
 *** Basic file I/O ***
 **********************/

static enum htsFormatCategory format_category(enum htsExactFormat fmt)
{
    switch (fmt) {
    case bam:
    case sam:
    case cram:
        return sequence_data;

    case vcf:
    case bcf:
        return variant_data;

    case bai:
    case crai:
    case csi:
    case gzi:
    case tbi:
        return index_file;

    case bed:
        return region_list;

    case htsget:
        return unknown_category;

    case unknown_format:
    case binary_format:
    case text_format:
    case format_maximum:
        break;
    }

    return unknown_category;
}

// Decompress up to ten or so bytes by peeking at the file, which must be
// positioned at the start of a GZIP block.
static size_t decompress_peek(hFILE *fp, unsigned char *dest, size_t destsize)
{
    // Typically at most a couple of hundred bytes of input are required
    // to get a few bytes of output from inflate(), so hopefully this buffer
    // size suffices in general.
    unsigned char buffer[512];
    z_stream zs;
    ssize_t npeek = hpeek(fp, buffer, sizeof buffer);

    if (npeek < 0) return 0;

    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = buffer;
    zs.avail_in = npeek;
    zs.next_out = dest;
    zs.avail_out = destsize;
    if (inflateInit2(&zs, 31) != Z_OK) return 0;

    while (zs.total_out < destsize)
        if (inflate(&zs, Z_SYNC_FLUSH) != Z_OK) break;

    destsize = zs.total_out;
    inflateEnd(&zs);

    return destsize;
}

// Parse "x.y" text, taking care because the string is not NUL-terminated
// and filling in major/minor only when the digits are followed by a delimiter,
// so we don't misread "1.10" as "1.1" due to reaching the end of the buffer.
static void
parse_version(htsFormat *fmt, const unsigned char *u, const unsigned char *ulim)
{
    const char *s    = (const char *) u;
    const char *slim = (const char *) ulim;
    short v;

    fmt->version.major = fmt->version.minor = -1;

    for (v = 0; s < slim && isdigit_c(*s); s++)
        v = 10 * v + *s - '0';

    if (s < slim) {
        fmt->version.major = v;
        if (*s == '.') {
            s++;
            for (v = 0; s < slim && isdigit_c(*s); s++)
                v = 10 * v + *s - '0';
            if (s < slim)
                fmt->version.minor = v;
        }
        else
            fmt->version.minor = 0;
    }
}

static int
cmp_nonblank(const char *key, const unsigned char *u, const unsigned char *ulim)
{
    const unsigned char *ukey = (const unsigned char *) key;

    while (*ukey)
        if (u >= ulim) return +1;
        else if (isspace_c(*u)) u++;
        else if (*u != *ukey) return (*ukey < *u)? -1 : +1;
        else u++, ukey++;

    return 0;
}

int hts_detect_format(hFILE *hfile, htsFormat *fmt)
{
    unsigned char s[32];
    ssize_t len = hpeek(hfile, s, 18);
    if (len < 0) return -1;

    if (len >= 2 && s[0] == 0x1f && s[1] == 0x8b) {
        // The stream is either gzip-compressed or BGZF-compressed.
        // Determine which, and decompress the first few bytes.
        fmt->compression = (len >= 18 && (s[3] & 4) &&
                            memcmp(&s[12], "BC\2\0", 4) == 0)? bgzf : gzip;
        len = decompress_peek(hfile, s, sizeof s);
    }
    else {
        fmt->compression = no_compression;
        len = hpeek(hfile, s, sizeof s);
    }
    if (len < 0) return -1;

    fmt->compression_level = -1;
    fmt->specific = NULL;

    if (len >= 6 && memcmp(s,"CRAM",4) == 0 && s[4]>=1 && s[4]<=3 && s[5]<=1) {
        fmt->category = sequence_data;
        fmt->format = cram;
        fmt->version.major = s[4], fmt->version.minor = s[5];
        fmt->compression = custom;
        return 0;
    }
    else if (len >= 4 && s[3] <= '\4') {
        if (memcmp(s, "BAM\1", 4) == 0) {
            fmt->category = sequence_data;
            fmt->format = bam;
            // TODO Decompress enough to pick version from @HD-VN header
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BAI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = bai;
            fmt->version.major = -1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BCF\4", 4) == 0) {
            fmt->category = variant_data;
            fmt->format = bcf;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BCF\2", 4) == 0) {
            fmt->category = variant_data;
            fmt->format = bcf;
            fmt->version.major = s[3];
            fmt->version.minor = (len >= 5 && s[4] <= 2)? s[4] : 0;
            return 0;
        }
        else if (memcmp(s, "CSI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = csi;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "TBI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = tbi;
            fmt->version.major = -1, fmt->version.minor = -1;
            return 0;
        }
    }
    else if (len >= 16 && memcmp(s, "##fileformat=VCF", 16) == 0) {
        fmt->category = variant_data;
        fmt->format = vcf;
        if (len >= 21 && s[16] == 'v')
            parse_version(fmt, &s[17], &s[len]);
        else
            fmt->version.major = fmt->version.minor = -1;
        return 0;
    }
    else if (len >= 4 && s[0] == '@' &&
             (memcmp(s, "@HD\t", 4) == 0 || memcmp(s, "@SQ\t", 4) == 0 ||
              memcmp(s, "@RG\t", 4) == 0 || memcmp(s, "@PG\t", 4) == 0)) {
        fmt->category = sequence_data;
        fmt->format = sam;
        // @HD-VN is not guaranteed to be the first tag, but then @HD is
        // not guaranteed to be present at all...
        if (len >= 9 && memcmp(s, "@HD\tVN:", 7) == 0)
            parse_version(fmt, &s[7], &s[len]);
        else
            fmt->version.major = 1, fmt->version.minor = -1;
        return 0;
    }
    else if (cmp_nonblank("{\"htsget\":", s, &s[len]) == 0) {
        fmt->category = unknown_category;
        fmt->format = htsget;
        fmt->version.major = fmt->version.minor = -1;
        return 0;
    }
    else {
        // Various possibilities for tab-delimited text:
        // .crai   (gzipped tab-delimited six columns: seqid 5*number)
        // .bed    ([3..12] tab-delimited columns)
        // .bedpe  (>= 10 tab-delimited columns)
        // .sam    (tab-delimited >= 11 columns: seqid number seqid...)
        // FIXME For now, assume it's SAM
        fmt->category = sequence_data;
        fmt->format = sam;
        fmt->version.major = 1, fmt->version.minor = -1;
        return 0;
    }

    fmt->category = unknown_category;
    fmt->format = unknown_format;
    fmt->version.major = fmt->version.minor = -1;
    fmt->compression = no_compression;
    return 0;
}

char *hts_format_description(const htsFormat *format)
{
    kstring_t str = { 0, 0, NULL };

    switch (format->format) {
    case sam:   kputs("SAM", &str); break;
    case bam:   kputs("BAM", &str); break;
    case cram:  kputs("CRAM", &str); break;
    case vcf:   kputs("VCF", &str); break;
    case bcf:
        if (format->version.major == 1) kputs("Legacy BCF", &str);
        else kputs("BCF", &str);
        break;
    case bai:   kputs("BAI", &str); break;
    case crai:  kputs("CRAI", &str); break;
    case csi:   kputs("CSI", &str); break;
    case tbi:   kputs("Tabix", &str); break;
    case htsget: kputs("htsget", &str); break;
    default:    kputs("unknown", &str); break;
    }

    if (format->version.major >= 0) {
        kputs(" version ", &str);
        kputw(format->version.major, &str);
        if (format->version.minor >= 0) {
            kputc('.', &str);
            kputw(format->version.minor, &str);
        }
    }

    switch (format->compression) {
    case custom: kputs(" compressed", &str); break;
    case gzip:   kputs(" gzip-compressed", &str); break;
    case bgzf:
        switch (format->format) {
        case bam:
        case bcf:
        case csi:
        case tbi:
            // These are by definition BGZF, so just use the generic term
            kputs(" compressed", &str);
            break;
        default:
            kputs(" BGZF-compressed", &str);
            break;
        }
        break;
    default: break;
    }

    switch (format->category) {
    case sequence_data: kputs(" sequence", &str); break;
    case variant_data:  kputs(" variant calling", &str); break;
    case index_file:    kputs(" index", &str); break;
    case region_list:   kputs(" genomic region", &str); break;
    default: break;
    }

    if (format->compression == no_compression)
        switch (format->format) {
        case sam:
        case crai:
        case vcf:
        case bed:
        case htsget:
            kputs(" text", &str);
            break;

        default:
            kputs(" data", &str);
            break;
        }
    else
        kputs(" data", &str);

    return ks_release(&str);
}

static htsFile *hts_open_format_impl(const char *fn, const char *mode, const htsFormat *fmt, hFILE* hf)
{
    char smode[102], *cp, *cp2, *mode_c;
    htsFile *fp = NULL;
    hFILE *hfile;
    char fmt_code = '\0';

    strncpy(smode, mode, 100);
    smode[100]=0;
    if ((cp = strchr(smode, ',')))
        *cp = '\0';

    // Migrate format code (b or c) to the end of the smode buffer.
    for (cp2 = cp = smode; *cp; cp++) {
        if (*cp == 'b')
            fmt_code = 'b';
        else if (*cp == 'c')
            fmt_code = 'c';
        else
            *cp2++ = *cp;
    }
    mode_c = cp2;
    *cp2++ = fmt_code;
    *cp2++ = 0;
    *cp2++ = 0;

    // Set or reset the format code if opts->format is used
    if (fmt && fmt->format != unknown_format)
        *mode_c = "\0g\0\0b\0c\0\0b\0g\0\0"[fmt->format];

    hfile = hf ? hf : hopen(fn, smode);
    if (hfile == NULL) goto error;

    fp = hts_hopen(hfile, fn, smode);
    if (fp == NULL) goto error;

    if (fmt && fmt->specific)
        if (hts_opt_apply(fp, fmt->specific) != 0)
            goto error;

    return fp;

error:
    hts_log_error("Failed to open file %s", fn);

    if (hfile)
        hclose_abruptly(hfile);

    return NULL;
}
htsFile *hts_open_format(const char *fn, const char *mode, const htsFormat *fmt)
{
	return hts_open_format_impl(fn, mode, fmt, NULL);
}

htsFile *hts_open_callback(const char* fn, hFILE_callback_ops* ops, const char* mode)
{
	if(NULL == ops) return NULL;
	hFILE* fp = hopen_callback(*ops, mode);
	return hts_open_format_impl(fn ? fn : "-", mode, NULL, fp);
}

htsFile *hts_open(const char *fn, const char *mode) {
    return hts_open_format(fn, mode, NULL);
}

/*
 * Splits str into a prefix, delimiter ('\0' or delim), and suffix, writing
 * the prefix in lowercase into buf and returning a pointer to the suffix.
 * On return, buf is always NUL-terminated; thus assumes that the "keyword"
 * prefix should be one of several known values of maximum length buflen-2.
 * (If delim is not found, returns a pointer to the '\0'.)
 */
static const char *
scan_keyword(const char *str, char delim, char *buf, size_t buflen)
{
    size_t i = 0;
    while (*str && *str != delim) {
        if (i < buflen-1) buf[i++] = tolower_c(*str);
        str++;
    }

    buf[i] = '\0';
    return *str? str+1 : str;
}

/*
 * Parses arg and appends it to the option list.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int hts_opt_add(hts_opt **opts, const char *c_arg) {
    hts_opt *o, *t;
    char *val;

    /*
     * IMPORTANT!!!
     * If you add another string option here, don't forget to also add
     * it to the case statement in hts_opt_apply.
     */

    if (!c_arg)
        return -1;

    if (!(o =  malloc(sizeof(*o))))
        return -1;

    if (!(o->arg = strdup(c_arg))) {
        free(o);
        return -1;
    }

    if (!(val = strchr(o->arg, '=')))
        val = "1"; // assume boolean
    else
        *val++ = '\0';

    if (strcmp(o->arg, "decode_md") == 0 ||
        strcmp(o->arg, "DECODE_MD") == 0)
        o->opt = CRAM_OPT_DECODE_MD, o->val.i = atoi(val);

    else if (strcmp(o->arg, "verbosity") == 0 ||
             strcmp(o->arg, "VERBOSITY") == 0)
        o->opt = CRAM_OPT_VERBOSITY, o->val.i = atoi(val);

    else if (strcmp(o->arg, "seqs_per_slice") == 0 ||
             strcmp(o->arg, "SEQS_PER_SLICE") == 0)
        o->opt = CRAM_OPT_SEQS_PER_SLICE, o->val.i = atoi(val);

    else if (strcmp(o->arg, "bases_per_slice") == 0 ||
             strcmp(o->arg, "BASES_PER_SLICE") == 0)
        o->opt = CRAM_OPT_BASES_PER_SLICE, o->val.i = atoi(val);

    else if (strcmp(o->arg, "slices_per_container") == 0 ||
             strcmp(o->arg, "SLICES_PER_CONTAINER") == 0)
        o->opt = CRAM_OPT_SLICES_PER_CONTAINER, o->val.i = atoi(val);

    else if (strcmp(o->arg, "embed_ref") == 0 ||
             strcmp(o->arg, "EMBED_REF") == 0)
        o->opt = CRAM_OPT_EMBED_REF, o->val.i = atoi(val);

    else if (strcmp(o->arg, "no_ref") == 0 ||
             strcmp(o->arg, "NO_REF") == 0)
        o->opt = CRAM_OPT_NO_REF, o->val.i = atoi(val);

    else if (strcmp(o->arg, "ignore_md5") == 0 ||
             strcmp(o->arg, "IGNORE_MD5") == 0)
        o->opt = CRAM_OPT_IGNORE_MD5, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_bzip2") == 0 ||
             strcmp(o->arg, "USE_BZIP2") == 0)
        o->opt = CRAM_OPT_USE_BZIP2, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_rans") == 0 ||
             strcmp(o->arg, "USE_RANS") == 0)
        o->opt = CRAM_OPT_USE_RANS, o->val.i = atoi(val);

    else if (strcmp(o->arg, "use_lzma") == 0 ||
             strcmp(o->arg, "USE_LZMA") == 0)
        o->opt = CRAM_OPT_USE_LZMA, o->val.i = atoi(val);

    else if (strcmp(o->arg, "reference") == 0 ||
             strcmp(o->arg, "REFERENCE") == 0)
        o->opt = CRAM_OPT_REFERENCE, o->val.s = val;

    else if (strcmp(o->arg, "version") == 0 ||
             strcmp(o->arg, "VERSION") == 0)
        o->opt = CRAM_OPT_VERSION, o->val.s =val;

    else if (strcmp(o->arg, "multi_seq_per_slice") == 0 ||
             strcmp(o->arg, "MULTI_SEQ_PER_SLICE") == 0)
        o->opt = CRAM_OPT_MULTI_SEQ_PER_SLICE, o->val.i = atoi(val);

    else if (strcmp(o->arg, "nthreads") == 0 ||
             strcmp(o->arg, "NTHREADS") == 0)
        o->opt = HTS_OPT_NTHREADS, o->val.i = atoi(val);

    else if (strcmp(o->arg, "cache_size") == 0 ||
             strcmp(o->arg, "CACHE_SIZE") == 0) {
        char *endp;
        o->opt = HTS_OPT_CACHE_SIZE;
        o->val.i = strtol(val, &endp, 0);
        // NB: Doesn't support floats, eg 1.5g
        // TODO: extend hts_parse_decimal? See also samtools sort.
        switch (*endp) {
        case 'g': case 'G': o->val.i *= 1024;
        case 'm': case 'M': o->val.i *= 1024;
        case 'k': case 'K': o->val.i *= 1024; break;
        case '\0': break;
        default:
            hts_log_error("Unrecognised cache size suffix '%c'", *endp);
            free(o->arg);
            free(o);
            return -1;
        }
    }

    else if (strcmp(o->arg, "required_fields") == 0 ||
             strcmp(o->arg, "REQUIRED_FIELDS") == 0)
        o->opt = CRAM_OPT_REQUIRED_FIELDS, o->val.i = strtol(val, NULL, 0);

    else if (strcmp(o->arg, "lossy_names") == 0 ||
             strcmp(o->arg, "LOSSY_NAMES") == 0)
        o->opt = CRAM_OPT_LOSSY_NAMES, o->val.i = strtol(val, NULL, 0);

    else if (strcmp(o->arg, "name_prefix") == 0 ||
             strcmp(o->arg, "NAME_PREFIX") == 0)
        o->opt = CRAM_OPT_PREFIX, o->val.s = val;

    else if (strcmp(o->arg, "store_md") == 0 ||
             strcmp(o->arg, "store_md") == 0)
        o->opt = CRAM_OPT_STORE_MD, o->val.i = atoi(val);

    else if (strcmp(o->arg, "store_nm") == 0 ||
             strcmp(o->arg, "store_nm") == 0)
        o->opt = CRAM_OPT_STORE_NM, o->val.i = atoi(val);

    else if (strcmp(o->arg, "block_size") == 0 ||
             strcmp(o->arg, "BLOCK_SIZE") == 0)
        o->opt = HTS_OPT_BLOCK_SIZE, o->val.i = strtol(val, NULL, 0);

    else if (strcmp(o->arg, "level") == 0 ||
             strcmp(o->arg, "LEVEL") == 0)
        o->opt = HTS_OPT_COMPRESSION_LEVEL, o->val.i = strtol(val, NULL, 0);

    else {
        hts_log_error("Unknown option '%s'", o->arg);
        free(o->arg);
        free(o);
        return -1;
    }

    o->next = NULL;

    // Append; assumes small list.
    if (*opts) {
        t = *opts;
        while (t->next)
            t = t->next;
        t->next = o;
    } else {
        *opts = o;
    }

    return 0;
}

/*
 * Applies an hts_opt option list to a given htsFile.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int hts_opt_apply(htsFile *fp, hts_opt *opts) {
    hts_opt *last = NULL;

    for (; opts;  opts = (last=opts)->next) {
        switch (opts->opt) {
            case CRAM_OPT_REFERENCE:
            case CRAM_OPT_VERSION:
            case CRAM_OPT_PREFIX:
                if (hts_set_opt(fp,  opts->opt,  opts->val.s) != 0)
                    return -1;
                break;
            default:
                if (hts_set_opt(fp,  opts->opt,  opts->val.i) != 0)
                    return -1;
                break;
        }
    }

    return 0;
}

/*
 * Frees an hts_opt list.
 */
void hts_opt_free(hts_opt *opts) {
    hts_opt *last = NULL;
    while (opts) {
        opts = (last=opts)->next;
        free(last->arg);
        free(last);
    }
}


/*
 * Tokenise options as (key(=value)?,)*(key(=value)?)?
 * NB: No provision for ',' appearing in the value!
 * Add backslashing rules?
 *
 * This could be used as part of a general command line option parser or
 * as a string concatenated onto the file open mode.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
int hts_parse_opt_list(htsFormat *fmt, const char *str) {
    while (str && *str) {
        const char *str_start;
        int len;
        char arg[8001];

        while (*str && *str == ',')
            str++;

        for (str_start = str; *str && *str != ','; str++);
        len = str - str_start;

        // Produce a nul terminated copy of the option
        strncpy(arg, str_start, len < 8000 ? len : 8000);
        arg[len < 8000 ? len : 8000] = '\0';

        if (hts_opt_add((hts_opt **)&fmt->specific, arg) != 0)
            return -1;

        if (*str)
            str++;
    }

    return 0;
}

/*
 * Accepts a string file format (sam, bam, cram, vcf, bam) optionally
 * followed by a comma separated list of key=value options and splits
 * these up into the fields of htsFormat struct.
 *
 * format is assumed to be already initialised, either to blank
 * "unknown" values or via previous hts_opt_add calls.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
int hts_parse_format(htsFormat *format, const char *str) {
    char fmt[8];
    const char *cp = scan_keyword(str, ',', fmt, sizeof fmt);

    format->version.minor = 0; // unknown
    format->version.major = 0; // unknown

    if (strcmp(fmt, "sam") == 0) {
        format->category          = sequence_data;
        format->format            = sam;
        format->compression       = no_compression;;
        format->compression_level = 0;
    } else if (strcmp(fmt, "bam") == 0) {
        format->category          = sequence_data;
        format->format            = bam;
        format->compression       = bgzf;
        format->compression_level = -1;
    } else if (strcmp(fmt, "cram") == 0) {
        format->category          = sequence_data;
        format->format            = cram;
        format->compression       = custom;
        format->compression_level = -1;
    } else if (strcmp(fmt, "vcf") == 0) {
        format->category          = variant_data;
        format->format            = vcf;
        format->compression       = no_compression;;
        format->compression_level = 0;
    } else if (strcmp(fmt, "bcf") == 0) {
        format->category          = variant_data;
        format->format            = bcf;
        format->compression       = bgzf;
        format->compression_level = -1;
    } else {
        return -1;
    }

    return hts_parse_opt_list(format, cp);
}


/*
 * Tokenise options as (key(=value)?,)*(key(=value)?)?
 * NB: No provision for ',' appearing in the value!
 * Add backslashing rules?
 *
 * This could be used as part of a general command line option parser or
 * as a string concatenated onto the file open mode.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
static int hts_process_opts(htsFile *fp, const char *opts) {
    htsFormat fmt;

    fmt.specific = NULL;
    if (hts_parse_opt_list(&fmt, opts) != 0)
        return -1;

    if (hts_opt_apply(fp, fmt.specific) != 0) {
        hts_opt_free(fmt.specific);
        return -1;
    }

    hts_opt_free(fmt.specific);

    return 0;
}


htsFile *hts_hopen(hFILE *hfile, const char *fn, const char *mode)
{
    hFILE *hfile_orig = hfile;
    htsFile *fp = (htsFile*)calloc(1, sizeof(htsFile));
    char simple_mode[101], *cp, *opts;
    simple_mode[100] = '\0';

    if (fp == NULL) goto error;

    fp->fn = strdup(fn);
    fp->is_be = ed_is_big();

    // Split mode into simple_mode,opts strings
    if ((cp = strchr(mode, ','))) {
        strncpy(simple_mode, mode, cp-mode <= 100 ? cp-mode : 100);
        simple_mode[cp-mode] = '\0';
        opts = cp+1;
    } else {
        strncpy(simple_mode, mode, 100);
        opts = NULL;
    }

    if (strchr(simple_mode, 'r')) {
        if (hts_detect_format(hfile, &fp->format) < 0) goto error;

        if (fp->format.format == htsget) {
            hFILE *hfile2 = hopen_htsget_redirect(hfile, simple_mode);
            if (hfile2 == NULL) goto error;

            // Build fp against the result of the redirection
            hfile = hfile2;
            if (hts_detect_format(hfile, &fp->format) < 0) goto error;
        }
    }
    else if (strchr(simple_mode, 'w') || strchr(simple_mode, 'a')) {
        htsFormat *fmt = &fp->format;
        fp->is_write = 1;

        if (strchr(simple_mode, 'b')) fmt->format = binary_format;
        else if (strchr(simple_mode, 'c')) fmt->format = cram;
        else fmt->format = text_format;

        if (strchr(simple_mode, 'z')) fmt->compression = bgzf;
        else if (strchr(simple_mode, 'g')) fmt->compression = gzip;
        else if (strchr(simple_mode, 'u')) fmt->compression = no_compression;
        else {
            // No compression mode specified, set to the default for the format
            switch (fmt->format) {
            case binary_format: fmt->compression = bgzf; break;
            case cram: fmt->compression = custom; break;
            case text_format: fmt->compression = no_compression; break;
            default: abort();
            }
        }

        // Fill in category (if determinable; e.g. 'b' could be BAM or BCF)
        fmt->category = format_category(fmt->format);

        fmt->version.major = fmt->version.minor = -1;
        fmt->compression_level = -1;
        fmt->specific = NULL;
    }
    else { errno = EINVAL; goto error; }

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        fp->fp.bgzf = bgzf_hopen(hfile, simple_mode);
        if (fp->fp.bgzf == NULL) goto error;
        fp->is_bin = fp->is_bgzf = 1;
        break;

    case cram:
        fp->fp.cram = cram_dopen(hfile, fn, simple_mode);
        if (fp->fp.cram == NULL) goto error;
        if (!fp->is_write)
            cram_set_option(fp->fp.cram, CRAM_OPT_DECODE_MD, 1);
        fp->is_cram = 1;
        break;

    case text_format:
    case sam:
    case vcf:
        if (fp->format.compression != no_compression) {
            fp->fp.bgzf = bgzf_hopen(hfile, simple_mode);
            if (fp->fp.bgzf == NULL) goto error;
            fp->is_bgzf = 1;
        }
        else
            fp->fp.hfile = hfile;
        break;

    default:
        errno = ENOEXEC;
        goto error;
    }

    if (opts)
        hts_process_opts(fp, opts);

    // If redirecting, close the original hFILE now (pedantically we would
    // instead close it in hts_close(), but this a simplifying optimisation)
    if (hfile != hfile_orig) hclose_abruptly(hfile_orig);

    return fp;

error:
    hts_log_error("Failed to open file %s", fn);

    // If redirecting, close the failed redirection hFILE that we have opened
    if (hfile != hfile_orig) hclose_abruptly(hfile);

    if (fp) {
        free(fp->fn);
        free(fp->fn_aux);
        free(fp);
    }
    return NULL;
}

int hts_close(htsFile *fp)
{
    int ret, save;

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        ret = bgzf_close(fp->fp.bgzf);
        break;

    case cram:
        if (!fp->is_write) {
            switch (cram_eof(fp->fp.cram)) {
            case 2:
                hts_log_warning("EOF marker is absent. The input is probably truncated");
                break;
            case 0:  /* not at EOF, but may not have wanted all seqs */
            default: /* case 1, expected EOF */
                break;
            }
        }
        ret = cram_close(fp->fp.cram);
        break;

    case text_format:
    case sam:
    case vcf:
        if (fp->format.compression != no_compression)
            ret = bgzf_close(fp->fp.bgzf);
        else
            ret = hclose(fp->fp.hfile);
        break;

    default:
        ret = -1;
        break;
    }

    save = errno;
    free(fp->fn);
    free(fp->fn_aux);
    free(fp->line.s);
    free(fp);
    errno = save;
    return ret;
}

const htsFormat *hts_get_format(htsFile *fp)
{
    return fp? &fp->format : NULL;
}

const char *hts_format_file_extension(const htsFormat *format) {
    if (!format)
        return "?";

    switch (format->format) {
    case sam:  return "sam";
    case bam:  return "bam";
    case bai:  return "bai";
    case cram: return "cram";
    case crai: return "crai";
    case vcf:  return "vcf";
    case bcf:  return "bcf";
    case csi:  return "csi";
    case gzi:  return "gzi";
    case tbi:  return "tbi";
    case bed:  return "bed";
    default:   return "?";
    }
}

static hFILE *hts_hfile(htsFile *fp) {
    switch (fp->format.format) {
    case binary_format: // fall through; still valid if bcf?
    case bam:          return bgzf_hfile(fp->fp.bgzf);
    case cram:         return cram_hfile(fp->fp.cram);
    case text_format:  return fp->fp.hfile;
    case sam:          return fp->fp.hfile;
    default:           return NULL;
    }
}

int hts_set_opt(htsFile *fp, enum hts_fmt_option opt, ...) {
    int r;
    va_list args;

    switch (opt) {
    case HTS_OPT_NTHREADS: {
        va_start(args, opt);
        int nthreads = va_arg(args, int);
        va_end(args);
        return hts_set_threads(fp, nthreads);
    }

    case HTS_OPT_BLOCK_SIZE: {
        hFILE *hf = hts_hfile(fp);

        if (hf) {
            va_start(args, opt);
            if (hfile_set_blksize(hf, va_arg(args, int)) != 0)
                hts_log_warning("Failed to change block size");
            va_end(args);
        }
        else {
            // To do - implement for vcf/bcf.
            hts_log_warning("Cannot change block size for this format");
        }

        return 0;
    }

    case HTS_OPT_THREAD_POOL: {
        va_start(args, opt);
        htsThreadPool *p = va_arg(args, htsThreadPool *);
        va_end(args);
        return hts_set_thread_pool(fp, p);
    }

    case HTS_OPT_CACHE_SIZE: {
        va_start(args, opt);
        int cache_size = va_arg(args, int);
        va_end(args);
        hts_set_cache_size(fp, cache_size);
        return 0;
    }

    case HTS_OPT_COMPRESSION_LEVEL: {
        va_start(args, opt);
        int level = va_arg(args, int);
        va_end(args);
        if (fp->is_bgzf)
            fp->fp.bgzf->compress_level = level;
    }

    default:
        break;
    }

    if (fp->format.format != cram)
        return 0;

    va_start(args, opt);
    r = cram_set_voption(fp->fp.cram, opt, args);
    va_end(args);

    return r;
}

BGZF *hts_get_bgzfp(htsFile *fp);

int hts_set_threads(htsFile *fp, int n)
{
    if (fp->format.compression == bgzf) {
        return bgzf_mt(hts_get_bgzfp(fp), n, 256/*unused*/);
    } else if (fp->format.format == cram) {
        return hts_set_opt(fp, CRAM_OPT_NTHREADS, n);
    }
    else return 0;
}

int hts_set_thread_pool(htsFile *fp, htsThreadPool *p) {
    if (fp->format.compression == bgzf) {
        return bgzf_thread_pool(hts_get_bgzfp(fp), p->pool, p->qsize);
    } else if (fp->format.format == cram) {
        return hts_set_opt(fp, CRAM_OPT_THREAD_POOL, p);
    }
    else return 0;
}

void hts_set_cache_size(htsFile *fp, int n)
{
    if (fp->format.compression == bgzf)
        bgzf_set_cache_size(hts_get_bgzfp(fp), n);
}

int hts_set_fai_filename(htsFile *fp, const char *fn_aux)
{
    free(fp->fn_aux);
    if (fn_aux) {
        fp->fn_aux = strdup(fn_aux);
        if (fp->fn_aux == NULL) return -1;
    }
    else fp->fn_aux = NULL;

    if (fp->format.format == cram)
        if (cram_set_option(fp->fp.cram, CRAM_OPT_REFERENCE, fp->fn_aux))
            return -1;

    return 0;
}

// For VCF/BCF backward sweeper. Not exposing these functions because their
// future is uncertain. Things will probably have to change with hFILE...
BGZF *hts_get_bgzfp(htsFile *fp)
{
    if (fp->is_bgzf)
        return fp->fp.bgzf;
    else
        return NULL;
}
int hts_useek(htsFile *fp, long uoffset, int where)
{
    if (fp->is_bgzf)
        return bgzf_useek(fp->fp.bgzf, uoffset, where);
    else
        return (hseek(fp->fp.hfile, uoffset, SEEK_SET) >= 0)? 0 : -1;
}
long hts_utell(htsFile *fp)
{
    if (fp->is_bgzf)
        return bgzf_utell(fp->fp.bgzf);
    else
        return htell(fp->fp.hfile);
}

int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
{
    int ret;
    if (! (delimiter == KS_SEP_LINE || delimiter == '\n')) {
        hts_log_error("Unexpected delimiter %d", delimiter);
        abort();
    }

    switch (fp->format.compression) {
    case no_compression:
        str->l = 0;
        ret = kgetline(str, (kgets_func *) hgets, fp->fp.hfile);
        if (ret >= 0) ret = str->l;
        else if (herrno(fp->fp.hfile)) ret = -2, errno = herrno(fp->fp.hfile);
        else ret = -1;
        break;

    case gzip:
    case bgzf:
        ret = bgzf_getline(fp->fp.bgzf, '\n', str);
        break;

    default:
        abort();
    }

    ++fp->lineno;
    return ret;
}

char **hts_readlist(const char *string, int is_file, int *_n)
{
    int m = 0, n = 0;
    char **s = 0;
    if ( is_file )
    {
        BGZF *fp = bgzf_open(string, "r");
        if ( !fp ) return NULL;

        kstring_t str;
        str.s = 0; str.l = str.m = 0;
        while (bgzf_getline(fp, '\n', &str) >= 0)
        {
            if (str.l == 0) continue;
            n++;
            hts_expand(char*,n,m,s);
            s[n-1] = strdup(str.s);
        }
        bgzf_close(fp);
        free(str.s);
    }
    else
    {
        const char *q = string, *p = string;
        while ( 1 )
        {
            if (*p == ',' || *p == 0)
            {
                n++;
                hts_expand(char*,n,m,s);
                s[n-1] = (char*)calloc(p - q + 1, 1);
                strncpy(s[n-1], q, p - q);
                q = p + 1;
            }
            if ( !*p ) break;
            p++;
        }
    }
    s = (char**)realloc(s, n * sizeof(char*));
    *_n = n;
    return s;
}

char **hts_readlines(const char *fn, int *_n)
{
    int m = 0, n = 0;
    char **s = 0;
    BGZF *fp = bgzf_open(fn, "r");
    if ( fp ) { // read from file
        kstring_t str;
        str.s = 0; str.l = str.m = 0;
        while (bgzf_getline(fp, '\n', &str) >= 0) {
            if (str.l == 0) continue;
            if (m == n) {
                m = m? m<<1 : 16;
                s = (char**)realloc(s, m * sizeof(char*));
            }
            s[n++] = strdup(str.s);
        }
        bgzf_close(fp);
        s = (char**)realloc(s, n * sizeof(char*));
        free(str.s);
    } else if (*fn == ':') { // read from string
        const char *q, *p;
        for (q = p = fn + 1;; ++p)
            if (*p == ',' || *p == 0) {
                if (m == n) {
                    m = m? m<<1 : 16;
                    s = (char**)realloc(s, m * sizeof(char*));
                }
                s[n] = (char*)calloc(p - q + 1, 1);
                strncpy(s[n++], q, p - q);
                q = p + 1;
                if (*p == 0) break;
            }
    } else return 0;
    s = (char**)realloc(s, n * sizeof(char*));
    *_n = n;
    return s;
}

// DEPRECATED: To be removed in a future HTSlib release
int hts_file_type(const char *fname)
{
    int len = strlen(fname);
    if ( !strcasecmp(".vcf.gz",fname+len-7) ) return FT_VCF_GZ;
    if ( !strcasecmp(".vcf",fname+len-4) ) return FT_VCF;
    if ( !strcasecmp(".bcf",fname+len-4) ) return FT_BCF_GZ;
    if ( !strcmp("-",fname) ) return FT_STDIN;

    hFILE *f = hopen(fname, "r");
    if (f == NULL) return 0;

    htsFormat fmt;
    if (hts_detect_format(f, &fmt) < 0) { hclose_abruptly(f); return 0; }
    if (hclose(f) < 0) return 0;

    switch (fmt.format) {
    case vcf: return (fmt.compression == no_compression)? FT_VCF : FT_VCF_GZ;
    case bcf: return (fmt.compression == no_compression)? FT_BCF : FT_BCF_GZ;
    default:  return 0;
    }
}

int hts_check_EOF(htsFile *fp)
{
    if (fp->format.compression == bgzf)
        return bgzf_check_EOF(hts_get_bgzfp(fp));
    else if (fp->format.format == cram)
        return cram_check_EOF(fp->fp.cram);
    else
        return 3;
}


/****************
 *** Indexing ***
 ****************/

#define HTS_MIN_MARKER_DIST 0x10000

// Finds the special meta bin
//  ((1<<(3 * n_lvls + 3)) - 1) / 7 + 1
#define META_BIN(idx) ((idx)->n_bins + 1)

#define pair64_lt(a,b) ((a).u < (b).u)

KSORT_INIT(_off, hts_pair64_t, pair64_lt)
KSORT_INIT(_off_max, hts_pair64_max_t, pair64_lt)

typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    int32_t n, m;
    uint64_t *offset;
} lidx_t;

struct __hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta; // MUST have a terminating NUL on the end
    struct {
        uint32_t last_bin, save_bin;
        int last_coor, last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};

static char * idx_format_name(int fmt) {
    switch (fmt) {
        case HTS_FMT_CSI: return "csi";
        case HTS_FMT_BAI: return "bai";
        case HTS_FMT_TBI: return "tbi";
        case HTS_FMT_CRAI: return "crai";
        default: return "unknown";
    }
}

static inline int insert_to_b(bidx_t *b, int bin, uint64_t beg, uint64_t end)
{
    khint_t k;
    bins_t *l;
    int absent;
    k = kh_put(bin, b, bin, &absent);
    if (absent < 0) return -1; // Out of memory
    l = &kh_value(b, k);
    if (absent) {
        l->m = 1; l->n = 0;
        l->list = (hts_pair64_t*)calloc(l->m, sizeof(hts_pair64_t));
        if (!l->list) {
            kh_del(bin, b, k);
            return -1;
        }
    } else if (l->n == l->m) {
        uint32_t new_m = l->m ? l->m << 1 : 1;
        hts_pair64_t *new_list = realloc(l->list, new_m * sizeof(hts_pair64_t));
        if (!new_list) return -1;
        l->list = new_list;
        l->m = new_m;
    }
    l->list[l->n].u = beg;
    l->list[l->n++].v = end;
    return 0;
}

static inline int insert_to_l(lidx_t *l, int64_t _beg, int64_t _end, uint64_t offset, int min_shift)
{
    int i, beg, end;
    beg = _beg >> min_shift;
    end = (_end - 1) >> min_shift;
    if (l->m < end + 1) {
        size_t new_m = l->m * 2 > end + 1 ? l->m * 2 : end + 1;
        uint64_t *new_offset;

        new_offset = (uint64_t*)realloc(l->offset, new_m * sizeof(uint64_t));
        if (!new_offset) return -1;

        // fill unused memory with (uint64_t)-1
        memset(new_offset + l->m, 0xff, sizeof(uint64_t) * (new_m - l->m));
        l->m = new_m;
        l->offset = new_offset;
    }
    for (i = beg; i <= end; ++i) {
        if (l->offset[i] == (uint64_t)-1) l->offset[i] = offset;
    }
    if (l->n < end + 1) l->n = end + 1;
    return 0;
}

hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls)
{
    hts_idx_t *idx;
    idx = (hts_idx_t*)calloc(1, sizeof(hts_idx_t));
    if (idx == NULL) return NULL;
    idx->fmt = fmt;
    idx->min_shift = min_shift;
    idx->n_lvls = n_lvls;
    idx->n_bins = ((1<<(3 * n_lvls + 3)) - 1) / 7;
    idx->z.save_bin = idx->z.save_tid = idx->z.last_tid = idx->z.last_bin = 0xffffffffu;
    idx->z.save_off = idx->z.last_off = idx->z.off_beg = idx->z.off_end = offset0;
    idx->z.last_coor = 0xffffffffu;
    if (n) {
        idx->n = idx->m = n;
        idx->bidx = (bidx_t**)calloc(n, sizeof(bidx_t*));
        if (idx->bidx == NULL) { free(idx); return NULL; }
        idx->lidx = (lidx_t*) calloc(n, sizeof(lidx_t));
        if (idx->lidx == NULL) { free(idx->bidx); free(idx); return NULL; }
    }
    return idx;
}

static void update_loff(hts_idx_t *idx, int i, int free_lidx)
{
    bidx_t *bidx = idx->bidx[i];
    lidx_t *lidx = &idx->lidx[i];
    khint_t k;
    int l;
    uint64_t offset0 = 0;
    if (bidx) {
        k = kh_get(bin, bidx, META_BIN(idx));
        if (k != kh_end(bidx))
            offset0 = kh_val(bidx, k).list[0].u;
        for (l = 0; l < lidx->n && lidx->offset[l] == (uint64_t)-1; ++l)
            lidx->offset[l] = offset0;
    } else l = 1;
    for (; l < lidx->n; ++l) // fill missing values
        if (lidx->offset[l] == (uint64_t)-1)
            lidx->offset[l] = lidx->offset[l-1];
    if (bidx == 0) return;
    for (k = kh_begin(bidx); k != kh_end(bidx); ++k) // set loff
        if (kh_exist(bidx, k))
        {
            if ( kh_key(bidx, k) < idx->n_bins )
            {
                int bot_bin = hts_bin_bot(kh_key(bidx, k), idx->n_lvls);
                // disable linear index if bot_bin out of bounds
                kh_val(bidx, k).loff = bot_bin < lidx->n ? lidx->offset[bot_bin] : 0;
            }
            else
                kh_val(bidx, k).loff = 0;
        }
    if (free_lidx) {
        free(lidx->offset);
        lidx->m = lidx->n = 0;
        lidx->offset = 0;
    }
}

static void compress_binning(hts_idx_t *idx, int i)
{
    bidx_t *bidx = idx->bidx[i];
    khint_t k;
    int l, m;
    if (bidx == 0) return;
    // merge a bin to its parent if the bin is too small
    for (l = idx->n_lvls; l > 0; --l) {
        unsigned start = hts_bin_first(l);
        for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
            bins_t *p, *q;
            if (!kh_exist(bidx, k) || kh_key(bidx, k) >= idx->n_bins || kh_key(bidx, k) < start) continue;
            p = &kh_value(bidx, k);
            if (l < idx->n_lvls && p->n > 1) ks_introsort(_off, p->n, p->list);
            if ((p->list[p->n - 1].v>>16) - (p->list[0].u>>16) < HTS_MIN_MARKER_DIST) {
                khint_t kp;
                kp = kh_get(bin, bidx, hts_bin_parent(kh_key(bidx, k)));
                if (kp == kh_end(bidx)) continue;
                q = &kh_val(bidx, kp);
                if (q->n + p->n > q->m) {
                    q->m = q->n + p->n;
                    kroundup32(q->m);
                    q->list = (hts_pair64_t*)realloc(q->list, q->m * sizeof(hts_pair64_t));
                }
                memcpy(q->list + q->n, p->list, p->n * sizeof(hts_pair64_t));
                q->n += p->n;
                free(p->list);
                kh_del(bin, bidx, k);
            }
        }
    }
    k = kh_get(bin, bidx, 0);
    if (k != kh_end(bidx)) ks_introsort(_off, kh_val(bidx, k).n, kh_val(bidx, k).list);
    // merge adjacent chunks that start from the same BGZF block
    for (k = kh_begin(bidx); k != kh_end(bidx); ++k) {
        bins_t *p;
        if (!kh_exist(bidx, k) || kh_key(bidx, k) >= idx->n_bins) continue;
        p = &kh_value(bidx, k);
        for (l = 1, m = 0; l < p->n; ++l) {
            if (p->list[m].v>>16 >= p->list[l].u>>16) {
                if (p->list[m].v < p->list[l].v) p->list[m].v = p->list[l].v;
            } else p->list[++m] = p->list[l];
        }
        p->n = m + 1;
    }
}

void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)
{
    int i;
    if (idx == NULL || idx->z.finished) return; // do not run this function on an empty index or multiple times
    if (idx->z.save_tid >= 0) {
        insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin, idx->z.save_off, final_offset);
        insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.off_beg, final_offset);
        insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx), idx->z.n_mapped, idx->z.n_unmapped);
    }
    for (i = 0; i < idx->n; ++i) {
        update_loff(idx, i, (idx->fmt == HTS_FMT_CSI));
        compress_binning(idx, i);
    }
    idx->z.finished = 1;
}

int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped)
{
    int bin;
    int64_t maxpos = (int64_t) 1 << (idx->min_shift + idx->n_lvls * 3);
    if (tid<0) beg = -1, end = 0;
    if (tid >= 0 && (beg > maxpos || end > maxpos)) {
        goto pos_too_big;
    }
    if (tid >= idx->m) { // enlarge the index
        uint32_t new_m = idx->m * 2 > tid + 1 ? idx->m * 2 : tid + 1;
        bidx_t **new_bidx;
        lidx_t *new_lidx;

        new_bidx = (bidx_t**)realloc(idx->bidx, new_m * sizeof(bidx_t*));
        if (!new_bidx) return -1;
        idx->bidx = new_bidx;

        new_lidx = (lidx_t*) realloc(idx->lidx, new_m * sizeof(lidx_t));
        if (!new_lidx) return -1;
        idx->lidx = new_lidx;

        memset(&idx->bidx[idx->m], 0, (new_m - idx->m) * sizeof(bidx_t*));
        memset(&idx->lidx[idx->m], 0, (new_m - idx->m) * sizeof(lidx_t));
        idx->m = new_m;
    }
    if (idx->n < tid + 1) idx->n = tid + 1;
    if (idx->z.finished) return 0;
    if (idx->z.last_tid != tid || (idx->z.last_tid >= 0 && tid < 0)) { // change of chromosome
        if ( tid>=0 && idx->n_no_coor )
        {
            hts_log_error("NO_COOR reads not in a single block at the end %d %d", tid, idx->z.last_tid);
            return -1;
        }
        if (tid>=0 && idx->bidx[tid] != 0)
        {
            hts_log_error("Chromosome blocks not continuous");
            return -1;
        }
        idx->z.last_tid = tid;
        idx->z.last_bin = 0xffffffffu;
    } else if (tid >= 0 && idx->z.last_coor > beg) { // test if positions are out of order
        hts_log_error("Unsorted positions on sequence #%d: %d followed by %d", tid+1, idx->z.last_coor+1, beg+1);
        return -1;
    }
    if ( tid>=0 )
    {
        if (idx->bidx[tid] == 0) idx->bidx[tid] = kh_init(bin);
        if (is_mapped) {
            // shoehorn [-1,0) (VCF POS=0) into the leftmost bottom-level bin
            if (beg < 0)  beg = 0;
            if (end <= 0) end = 1;
            // idx->z.last_off points to the start of the current record
            if (insert_to_l(&idx->lidx[tid], beg, end,
                            idx->z.last_off, idx->min_shift) < 0) return -1;
        }
    }
    else idx->n_no_coor++;
    bin = hts_reg2bin(beg, end, idx->min_shift, idx->n_lvls);
    if ((int)idx->z.last_bin != bin) { // then possibly write the binning index
        if (idx->z.save_bin != 0xffffffffu) { // save_bin==0xffffffffu only happens to the first record
            if (insert_to_b(idx->bidx[idx->z.save_tid], idx->z.save_bin,
                            idx->z.save_off, idx->z.last_off) < 0) return -1;
        }
        if (idx->z.last_bin == 0xffffffffu && idx->z.save_bin != 0xffffffffu) { // change of chr; keep meta information
            idx->z.off_end = idx->z.last_off;
            if (insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx),
                            idx->z.off_beg, idx->z.off_end) < 0) return -1;
            if (insert_to_b(idx->bidx[idx->z.save_tid], META_BIN(idx),
                            idx->z.n_mapped, idx->z.n_unmapped) < 0) return -1;
            idx->z.n_mapped = idx->z.n_unmapped = 0;
            idx->z.off_beg = idx->z.off_end;
        }
        idx->z.save_off = idx->z.last_off;
        idx->z.save_bin = idx->z.last_bin = bin;
        idx->z.save_tid = tid;
    }
    if (is_mapped) ++idx->z.n_mapped;
    else ++idx->z.n_unmapped;
    idx->z.last_off = offset;
    idx->z.last_coor = beg;
    return 0;

 pos_too_big: {
        int64_t max = end > beg ? end : beg, s = 1 << 14;
        int n_lvls = 0;
        while (max > s) {
            n_lvls++;
            s <<= 3;
        }

        if (idx->fmt == HTS_FMT_CSI) {
            hts_log_error("Region %d..%d cannot be stored in a csi index "
                "with min_shift = %d, n_lvls = %d. Try using "
                "min_shift = 14, n_lvls >= %d",
                beg, end,
                idx->min_shift, idx->n_lvls,
                n_lvls);
        } else {
            hts_log_error("Region %d..%d cannot be stored in a %s index. "
                "Try using a csi index with min_shift = 14, "
                "n_lvls >= %d",
                beg, end, idx_format_name(idx->fmt),
                n_lvls);
        }
        errno = ERANGE;
        return -1;
    }
}

void hts_idx_destroy(hts_idx_t *idx)
{
    khint_t k;
    int i;
    if (idx == 0) return;

    // For HTS_FMT_CRAI, idx actually points to a different type -- see sam.c
    if (idx->fmt == HTS_FMT_CRAI) {
        hts_cram_idx_t *cidx = (hts_cram_idx_t *) idx;
        cram_index_free(cidx->cram);
        free(cidx);
        return;
    }

    for (i = 0; i < idx->m; ++i) {
        bidx_t *bidx = idx->bidx[i];
        free(idx->lidx[i].offset);
        if (bidx == 0) continue;
        for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
            if (kh_exist(bidx, k))
                free(kh_value(bidx, k).list);
        kh_destroy(bin, bidx);
    }
    free(idx->bidx); free(idx->lidx); free(idx->meta);
    free(idx);
}

// The optimizer eliminates these ed_is_big() calls; still it would be good to
// TODO Determine endianness at configure- or compile-time

static inline ssize_t HTS_RESULT_USED idx_write_int32(BGZF *fp, int32_t x)
{
    if (ed_is_big()) x = ed_swap_4(x);
    return bgzf_write(fp, &x, sizeof x);
}

static inline ssize_t HTS_RESULT_USED idx_write_uint32(BGZF *fp, uint32_t x)
{
    if (ed_is_big()) x = ed_swap_4(x);
    return bgzf_write(fp, &x, sizeof x);
}

static inline ssize_t HTS_RESULT_USED idx_write_uint64(BGZF *fp, uint64_t x)
{
    if (ed_is_big()) x = ed_swap_8(x);
    return bgzf_write(fp, &x, sizeof x);
}

static inline void swap_bins(bins_t *p)
{
    int i;
    for (i = 0; i < p->n; ++i) {
        ed_swap_8p(&p->list[i].u);
        ed_swap_8p(&p->list[i].v);
    }
}

static int hts_idx_save_core(const hts_idx_t *idx, BGZF *fp, int fmt)
{
    int32_t i, j;

    #define check(ret) if ((ret) < 0) return -1

    check(idx_write_int32(fp, idx->n));
    if (fmt == HTS_FMT_TBI && idx->l_meta)
        check(bgzf_write(fp, idx->meta, idx->l_meta));

    for (i = 0; i < idx->n; ++i) {
        khint_t k;
        bidx_t *bidx = idx->bidx[i];
        lidx_t *lidx = &idx->lidx[i];
        // write binning index
        check(idx_write_int32(fp, bidx? kh_size(bidx) : 0));
        if (bidx)
            for (k = kh_begin(bidx); k != kh_end(bidx); ++k)
                if (kh_exist(bidx, k)) {
                    bins_t *p = &kh_value(bidx, k);
                    check(idx_write_uint32(fp, kh_key(bidx, k)));
                    if (fmt == HTS_FMT_CSI) check(idx_write_uint64(fp, p->loff));
                    //int j;for(j=0;j<p->n;++j)fprintf(stderr,"%d,%llx,%d,%llx:%llx\n",kh_key(bidx,k),kh_val(bidx, k).loff,j,p->list[j].u,p->list[j].v);
                    check(idx_write_int32(fp, p->n));
                    for (j = 0; j < p->n; ++j) {
                        check(idx_write_uint64(fp, p->list[j].u));
                        check(idx_write_uint64(fp, p->list[j].v));
                    }
                }

        // write linear index
        if (fmt != HTS_FMT_CSI) {
            check(idx_write_int32(fp, lidx->n));
            for (j = 0; j < lidx->n; ++j)
                check(idx_write_uint64(fp, lidx->offset[j]));
        }
    }

    check(idx_write_uint64(fp, idx->n_no_coor));
    return 0;
    #undef check
}

int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)
{
    int ret, save;
    char *fnidx = (char*)calloc(1, strlen(fn) + 5);
    if (fnidx == NULL) return -1;

    strcpy(fnidx, fn);
    switch (fmt) {
    case HTS_FMT_BAI: strcat(fnidx, ".bai"); break;
    case HTS_FMT_CSI: strcat(fnidx, ".csi"); break;
    case HTS_FMT_TBI: strcat(fnidx, ".tbi"); break;
    default: abort();
    }

    ret = hts_idx_save_as(idx, fn, fnidx, fmt);
    save = errno;
    free(fnidx);
    errno = save;
    return ret;
}

int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt)
{
    BGZF *fp;

    #define check(ret) if ((ret) < 0) goto fail

    if (fnidx == NULL) return hts_idx_save(idx, fn, fmt);

    fp = bgzf_open(fnidx, (fmt == HTS_FMT_BAI)? "wu" : "w");
    if (fp == NULL) return -1;

    if (fmt == HTS_FMT_CSI) {
        check(bgzf_write(fp, "CSI\1", 4));
        check(idx_write_int32(fp, idx->min_shift));
        check(idx_write_int32(fp, idx->n_lvls));
        check(idx_write_uint32(fp, idx->l_meta));
        if (idx->l_meta) check(bgzf_write(fp, idx->meta, idx->l_meta));
    } else if (fmt == HTS_FMT_TBI) {
        check(bgzf_write(fp, "TBI\1", 4));
    } else if (fmt == HTS_FMT_BAI) {
        check(bgzf_write(fp, "BAI\1", 4));
    } else abort();

    check(hts_idx_save_core(idx, fp, fmt));

    return bgzf_close(fp);
    #undef check

fail:
    bgzf_close(fp);
    return -1;
}

static int hts_idx_load_core(hts_idx_t *idx, BGZF *fp, int fmt)
{
    int32_t i, n, is_be;
    is_be = ed_is_big();
    if (idx == NULL) return -4;
    for (i = 0; i < idx->n; ++i) {
        bidx_t *h;
        lidx_t *l = &idx->lidx[i];
        uint32_t key;
        int j, absent;
        bins_t *p;
        h = idx->bidx[i] = kh_init(bin);
        if (bgzf_read(fp, &n, 4) != 4) return -1;
        if (is_be) ed_swap_4p(&n);
        for (j = 0; j < n; ++j) {
            khint_t k;
            if (bgzf_read(fp, &key, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&key);
            k = kh_put(bin, h, key, &absent);
            if (absent <= 0) return -3; // Duplicate bin number
            p = &kh_val(h, k);
            if (fmt == HTS_FMT_CSI) {
                if (bgzf_read(fp, &p->loff, 8) != 8) return -1;
                if (is_be) ed_swap_8p(&p->loff);
            } else p->loff = 0;
            if (bgzf_read(fp, &p->n, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&p->n);
            p->m = p->n;
            p->list = (hts_pair64_t*)malloc(p->m * sizeof(hts_pair64_t));
            if (p->list == NULL) return -2;
            if (bgzf_read(fp, p->list, p->n<<4) != p->n<<4) return -1;
            if (is_be) swap_bins(p);
        }
        if (fmt != HTS_FMT_CSI) { // load linear index
            int j;
            if (bgzf_read(fp, &l->n, 4) != 4) return -1;
            if (is_be) ed_swap_4p(&l->n);
            l->m = l->n;
            l->offset = (uint64_t*)malloc(l->n * sizeof(uint64_t));
            if (l->offset == NULL) return -2;
            if (bgzf_read(fp, l->offset, l->n << 3) != l->n << 3) return -1;
            if (is_be) for (j = 0; j < l->n; ++j) ed_swap_8p(&l->offset[j]);
            for (j = 1; j < l->n; ++j) // fill missing values; may happen given older samtools and tabix
                if (l->offset[j] == 0) l->offset[j] = l->offset[j-1];
            update_loff(idx, i, 1);
        }
    }
    if (bgzf_read(fp, &idx->n_no_coor, 8) != 8) idx->n_no_coor = 0;
    if (is_be) ed_swap_8p(&idx->n_no_coor);
    return 0;
}

static hts_idx_t *hts_idx_load_local(const char *fn)
{
    uint8_t magic[4];
    int i, is_be;
    hts_idx_t *idx = NULL;
    uint8_t *meta = NULL;
    BGZF *fp = bgzf_open(fn, "r");
    if (fp == NULL) return NULL;
    is_be = ed_is_big();
    if (bgzf_read(fp, magic, 4) != 4) goto fail;

    if (memcmp(magic, "CSI\1", 4) == 0) {
        uint32_t x[3], n;
        if (bgzf_read(fp, x, 12) != 12) goto fail;
        if (is_be) for (i = 0; i < 3; ++i) ed_swap_4p(&x[i]);
        if (x[2]) {
            if (SIZE_MAX - x[2] < 1) goto fail; // Prevent possible overflow
            if ((meta = (uint8_t*)malloc((size_t) x[2] + 1)) == NULL) goto fail;
            if (bgzf_read(fp, meta, x[2]) != x[2]) goto fail;
            // Prevent possible strlen past the end in tbx_index_load2
            meta[x[2]] = '\0';
        }
        if (bgzf_read(fp, &n, 4) != 4) goto fail;
        if (is_be) ed_swap_4p(&n);
        if ((idx = hts_idx_init(n, HTS_FMT_CSI, 0, x[0], x[1])) == NULL) goto fail;
        idx->l_meta = x[2];
        idx->meta = meta;
        meta = NULL;
        if (hts_idx_load_core(idx, fp, HTS_FMT_CSI) < 0) goto fail;
    }
    else if (memcmp(magic, "TBI\1", 4) == 0) {
        uint8_t x[8 * 4];
        uint32_t n;
        // Read file header
        if (bgzf_read(fp, x, sizeof(x)) != sizeof(x)) goto fail;
        n = le_to_u32(&x[0]); // location of n_ref
        if ((idx = hts_idx_init(n, HTS_FMT_TBI, 0, 14, 5)) == NULL) goto fail;
        n = le_to_u32(&x[7*4]); // location of l_nm
        if (n > UINT32_MAX - 29) goto fail; // Prevent possible overflow
        idx->l_meta = 28 + n;
        if ((idx->meta = (uint8_t*)malloc(idx->l_meta + 1)) == NULL) goto fail;
        // copy format, col_seq, col_beg, col_end, meta, skip, l_nm
        // N.B. left in little-endian byte order.
        memcpy(idx->meta, &x[1*4], 28);
        // Read in sequence names.
        if (bgzf_read(fp, idx->meta + 28, n) != n) goto fail;
        // Prevent possible strlen past the end in tbx_index_load2
        idx->meta[idx->l_meta] = '\0';
        if (hts_idx_load_core(idx, fp, HTS_FMT_TBI) < 0) goto fail;
    }
    else if (memcmp(magic, "BAI\1", 4) == 0) {
        uint32_t n;
        if (bgzf_read(fp, &n, 4) != 4) goto fail;
        if (is_be) ed_swap_4p(&n);
        idx = hts_idx_init(n, HTS_FMT_BAI, 0, 14, 5);
        if (hts_idx_load_core(idx, fp, HTS_FMT_BAI) < 0) goto fail;
    }
    else { errno = EINVAL; goto fail; }

    bgzf_close(fp);
    return idx;

fail:
    bgzf_close(fp);
    hts_idx_destroy(idx);
    free(meta);
    return NULL;
}

int hts_idx_set_meta(hts_idx_t *idx, uint32_t l_meta, uint8_t *meta,
                      int is_copy)
{
    uint8_t *new_meta = meta;
    if (is_copy) {
        size_t l = l_meta;
        if (l > SIZE_MAX - 1) {
            errno = ENOMEM;
            return -1;
        }
        new_meta = malloc(l + 1);
        if (!new_meta) return -1;
        memcpy(new_meta, meta, l);
        // Prevent possible strlen past the end in tbx_index_load2
        meta[l + 1] = '\0';
    }
    if (idx->meta) free(idx->meta);
    idx->l_meta = l_meta;
    idx->meta = new_meta;
    return 0;
}

uint8_t *hts_idx_get_meta(hts_idx_t *idx, uint32_t *l_meta)
{
    *l_meta = idx->l_meta;
    return idx->meta;
}

const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr)
{
    if ( !idx->n )
    {
        *n = 0;
        return NULL;
    }

    int tid = 0, i;
    const char **names = (const char**) calloc(idx->n,sizeof(const char*));
    for (i=0; i<idx->n; i++)
    {
        bidx_t *bidx = idx->bidx[i];
        if ( !bidx ) continue;
        names[tid++] = getid(hdr,i);
    }
    *n = tid;
    return names;
}

int hts_idx_get_stat(const hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped)
{
    if ( idx->fmt == HTS_FMT_CRAI ) {
        *mapped = 0; *unmapped = 0;
        return -1;
    }

    bidx_t *h = idx->bidx[tid];
    khint_t k = kh_get(bin, h, META_BIN(idx));
    if (k != kh_end(h)) {
        *mapped = kh_val(h, k).list[1].u;
        *unmapped = kh_val(h, k).list[1].v;
        return 0;
    } else {
        *mapped = 0; *unmapped = 0;
        return -1;
    }
}

uint64_t hts_idx_get_n_no_coor(const hts_idx_t* idx)
{
    return idx->n_no_coor;
}

/****************
 *** Iterator ***
 ****************/

static inline int reg2bins(int64_t beg, int64_t end, hts_itr_t *itr, int min_shift, int n_lvls)
{
    int l, t, s = min_shift + (n_lvls<<1) + n_lvls;
    if (beg >= end) return 0;
    if (end >= 1LL<<s) end = 1LL<<s;
    for (--end, l = 0, t = 0; l <= n_lvls; s -= 3, t += 1<<((l<<1)+l), ++l) {
        int b, e, n, i;
        b = t + (beg>>s); e = t + (end>>s); n = e - b + 1;
        if (itr->bins.n + n > itr->bins.m) {
            itr->bins.m = itr->bins.n + n;
            kroundup32(itr->bins.m);
            itr->bins.a = (int*)realloc(itr->bins.a, sizeof(int) * itr->bins.m);
        }
        for (i = b; i <= e; ++i) itr->bins.a[itr->bins.n++] = i;
    }
    return itr->bins.n;
}

static inline int reg2intervals(hts_itr_multi_t *iter, const hts_idx_t *idx, int tid, int64_t beg, int64_t end, uint64_t min_off, uint64_t max_off, int min_shift, int n_lvls)
{
    int l, t, s;
    int b, e, i, j;
    hts_pair64_max_t *off;
    bidx_t *bidx;
    khint_t k;

    if (!iter || !idx || (bidx = idx->bidx[tid]) == NULL || beg >= end)
        return -1;

    s = min_shift + (n_lvls<<1) + n_lvls;
    if (end >= 1LL<<s)
        end = 1LL<<s;

    for (--end, l = 0, t = 0; l <= n_lvls; s -= 3, t += 1<<((l<<1)+l), ++l) {
        b = t + (beg>>s); e = t + (end>>s);

        for (i = b; i <= e; ++i) {
            if ((k = kh_get(bin, bidx, i)) != kh_end(bidx)) {
                bins_t *p = &kh_value(bidx, k);

                if (p->n) {
                    off = (hts_pair64_max_t*)realloc(iter->off, (iter->n_off + p->n) * sizeof(hts_pair64_max_t));
                    if (!off)
                        return -2;

                    iter->off = off;
                    for (j = 0; j < p->n; ++j) {
                        if (p->list[j].v > min_off && p->list[j].u < max_off) {
                            iter->off[iter->n_off].u = p->list[j].u;
                            iter->off[iter->n_off].v = p->list[j].v;
                            iter->off[iter->n_off].max = ((uint64_t)tid<<32) | (end+1);
                            iter->n_off++;
                        }
                    }
                }
            }
        }
    }

    return iter->n_off;
}

static int compare_regions(const void *r1, const void *r2) {
    hts_reglist_t *reg1 = (hts_reglist_t *)r1;
    hts_reglist_t *reg2 = (hts_reglist_t *)r2;

    if (reg1->tid < 0 && reg2->tid >= 0)
        return 1;
    else if (reg1->tid >= 0 && reg2->tid < 0)
        return -1;
    else
        return reg1->tid - reg2->tid;
}

uint64_t hts_itr_off(const hts_idx_t* idx, int tid) {

    int i;
    bidx_t* bidx;
    uint64_t off0 = (uint64_t) -1;
    khint_t k;
    switch (tid) {
    case HTS_IDX_START:
        // Find the smallest offset, note that sequence ids may not be ordered sequentially
        for (i = 0; i < idx->n; i++) {
            bidx = idx->bidx[i];
            k = kh_get(bin, bidx, META_BIN(idx));
            if (k == kh_end(bidx))
                continue;

            if (off0 > kh_val(bidx, k).list[0].u)
                off0 = kh_val(bidx, k).list[0].u;
        }
        if (off0 == (uint64_t) -1 && idx->n_no_coor)
            off0 = 0;
        // only no-coor reads in this bam
        break;
    case HTS_IDX_NOCOOR:
        /* No-coor reads sort after all of the mapped reads.  The position
           is not stored in the index itself, so need to find the end
           offset for the last mapped read.  A loop is needed here in
           case references at the end of the file have no mapped reads,
           or sequence ids are not ordered sequentially.
           See issue samtools#568 and commits b2aab8, 60c22d and cc207d. */
        for (i = 0; i < idx->n; i++) {
            bidx = idx->bidx[i];
            k = kh_get(bin, bidx, META_BIN(idx));
            if (k != kh_end(bidx)) {
                if (off0 == (uint64_t) -1 || off0 < kh_val(bidx, k).list[0].v) {
                    off0 = kh_val(bidx, k).list[0].v;
                }
            }
        }
        if (off0 == (uint64_t) -1 && idx->n_no_coor)
            off0 = 0;
        // only no-coor reads in this bam
        break;
    case HTS_IDX_REST:
        off0 = 0;
        break;
    case HTS_IDX_NONE:
        off0 = 0;
        break;
    }

    return off0;
}

hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec)
{
    int i, n_off, l, bin;
    hts_pair64_t *off;
    khint_t k;
    bidx_t *bidx;
    uint64_t min_off, max_off;
    hts_itr_t *iter = (hts_itr_t*)calloc(1, sizeof(hts_itr_t));
    if (iter) {
        if (tid < 0) {
            uint64_t off = hts_itr_off(idx, tid);
            if (off != (uint64_t) -1) {
                iter->read_rest = 1;
                iter->curr_off = off;
                iter->readrec = readrec;
                if (tid == HTS_IDX_NONE)
                    iter->finished = 1;
            } else {
                free(iter);
                iter = NULL;
            }
        } else {
            if (beg < 0) beg = 0;
            if (end < beg) return 0;
            if (tid >= idx->n || (bidx = idx->bidx[tid]) == NULL) return 0;

            iter->tid = tid, iter->beg = beg, iter->end = end; iter->i = -1;
            iter->readrec = readrec;

            if ( !kh_size(bidx) ) { iter->finished = 1; return iter; }

            // compute min_off
            bin = hts_bin_first(idx->n_lvls) + (beg>>idx->min_shift);
            do {
                int first;
                k = kh_get(bin, bidx, bin);
                if (k != kh_end(bidx)) break;
                first = (hts_bin_parent(bin)<<3) + 1;
                if (bin > first) --bin;
                else bin = hts_bin_parent(bin);
            } while (bin);
            if (bin == 0) k = kh_get(bin, bidx, bin);
            min_off = k != kh_end(bidx)? kh_val(bidx, k).loff : 0;

            // compute max_off: a virtual offset from a bin to the right of end
            bin = hts_bin_first(idx->n_lvls) + ((end-1) >> idx->min_shift) + 1;
            if (bin >= idx->n_bins) bin = 0;
            while (1) {
                // search for an extant bin by moving right, but moving up to the
                // parent whenever we get to a first child (which also covers falling
                // off the RHS, which wraps around and immediately goes up to bin 0)
                while (bin % 8 == 1) bin = hts_bin_parent(bin);
                if (bin == 0) { max_off = (uint64_t)-1; break; }
                k = kh_get(bin, bidx, bin);
                if (k != kh_end(bidx) && kh_val(bidx, k).n > 0) { max_off = kh_val(bidx, k).list[0].u; break; }
                bin++;
            }

            // retrieve bins
            reg2bins(beg, end, iter, idx->min_shift, idx->n_lvls);

            for (i = n_off = 0; i < iter->bins.n; ++i)
                if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx))
                    n_off += kh_value(bidx, k).n;
            if (n_off == 0) {
                // No overlapping bins means the iterator has already finished.
                iter->finished = 1;
                return iter;
            }
            off = (hts_pair64_t*)calloc(n_off, sizeof(hts_pair64_t));
            for (i = n_off = 0; i < iter->bins.n; ++i) {
                if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx)) {
                    int j;
                    bins_t *p = &kh_value(bidx, k);
                    for (j = 0; j < p->n; ++j)
                        if (p->list[j].v > min_off && p->list[j].u < max_off)
                            off[n_off++] = p->list[j];
                }
            }

            if (n_off == 0) {
                free(off);
                iter->finished = 1;
                return iter;
            }
            ks_introsort(_off, n_off, off);
            // resolve completely contained adjacent blocks
            for (i = 1, l = 0; i < n_off; ++i)
                if (off[l].v < off[i].v) off[++l] = off[i];
            n_off = l + 1;
            // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
            for (i = 1; i < n_off; ++i)
                if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
            // merge adjacent blocks
            for (i = 1, l = 0; i < n_off; ++i) {
                if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
                else off[++l] = off[i];
            }
            n_off = l + 1;
            iter->n_off = n_off; iter->off = off;
        }
    }
    return iter;
}

hts_itr_multi_t *hts_itr_multi_bam(const hts_idx_t *idx, hts_itr_multi_t *iter)
{
    int i, j, l, n_off = 0, bin;
    hts_pair64_max_t *off = NULL;
    khint_t k;
    bidx_t *bidx;
    uint64_t min_off, max_off, t_off = (uint64_t)-1;
    int tid, beg, end;
    hts_reglist_t *curr_reg;

    if (iter) {
        iter->i = -1;
        for (i=0; i<iter->n_reg; i++) {

            curr_reg = &iter->reg_list[i];
            tid = curr_reg->tid;

            if (tid < 0) {
                t_off = hts_itr_off(idx, tid);
                if (t_off != (uint64_t)-1) {
                    switch (tid) {
                        case HTS_IDX_NONE:
                            iter->finished = 1;
                        case HTS_IDX_START:
                        case HTS_IDX_REST:
                            iter->curr_off = t_off;
                            iter->n_reg = 0;
                            iter->reg_list = NULL;
                            iter->read_rest = 1;
                            return iter;
                        case HTS_IDX_NOCOOR:
                            iter->nocoor = 1;
                            iter->nocoor_off = t_off;
                    }
                }
            } else {
                if (tid >= idx->n || (bidx = idx->bidx[tid]) == NULL || !kh_size(bidx))
                    continue;

                for(j=0; j<curr_reg->count; j++) {
                    hts_pair32_t *curr_intv = &curr_reg->intervals[j];
                    if (curr_intv->end < curr_intv->beg)
                        continue;

                    beg = curr_intv->beg;
                    end = curr_intv->end;

                    /* Compute 'min_off' by searching the lowest level bin containing 'beg'.
                       If the computed bin is not in the index, try the next bin to the
                       left, belonging to the same parent. If it is the first sibling bin,
                       try the parent bin. */
                    bin = hts_bin_first(idx->n_lvls) + (beg>>idx->min_shift);
                    do {
                        int first;
                        k = kh_get(bin, bidx, bin);
                        if (k != kh_end(bidx)) break;
                        first = (hts_bin_parent(bin)<<3) + 1;
                        if (bin > first) --bin;
                        else bin = hts_bin_parent(bin);
                    } while (bin);
                    if (bin == 0)
                        k = kh_get(bin, bidx, bin);
                    min_off = k != kh_end(bidx)? kh_val(bidx, k).loff : 0;

                    // compute max_off: a virtual offset from a bin to the right of end
                    bin = hts_bin_first(idx->n_lvls) + ((end-1) >> idx->min_shift) + 1;
                    if (bin >= idx->n_bins) bin = 0;
                    while (1) {
                    // search for an extant bin by moving right, but moving up to the
                    // parent whenever we get to a first child (which also covers falling
                    // off the RHS, which wraps around and immediately goes up to bin 0)
                        while (bin % 8 == 1) bin = hts_bin_parent(bin);
                        if (bin == 0) { max_off = (uint64_t)-1; break; }
                        k = kh_get(bin, bidx, bin);
                        if (k != kh_end(bidx) && kh_val(bidx, k).n > 0) {
                            max_off = kh_val(bidx, k).list[0].u;
                            break;
                        }
                        bin++;
                    }

                    //convert coordinates to file offsets
                    reg2intervals(iter, idx, tid, beg, end, min_off, max_off, idx->min_shift, idx->n_lvls);
                }
            }
        }

        off = iter->off;
        n_off = iter->n_off;

        if (n_off) {
            ks_introsort(_off_max, n_off, off);
            // resolve completely contained adjacent blocks
            for (i = 1, l = 0; i < n_off; ++i) {
                if (off[l].v < off[i].v) {
                    off[++l] = off[i];
                } else {
                    off[l].max = (off[i].max > off[l].max ? off[i].max : off[l].max);
                }
            }
            n_off = l + 1;
            // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
            for (i = 1; i < n_off; ++i)
                if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
            // merge adjacent blocks
            for (i = 1, l = 0; i < n_off; ++i) {
                if (off[l].v>>16 == off[i].u>>16) {
                    off[l].v = off[i].v;
                    off[l].max = (off[i].max > off[l].max ? off[i].max : off[l].max);
                } else off[++l] = off[i];
            }
            n_off = l + 1;
            iter->n_off = n_off; iter->off = off;
        }

        if(!n_off && !iter->nocoor)
            iter->finished = 1;
    }
    return iter;
}

hts_itr_multi_t *hts_itr_multi_cram(const hts_idx_t *idx, hts_itr_multi_t *iter)
{
    const hts_cram_idx_t *cidx = (const hts_cram_idx_t *) idx;
    int tid, beg, end, i, j, l, n_off = 0;
    hts_reglist_t *curr_reg;
    hts_pair32_t *curr_intv;
    hts_pair64_max_t *off = NULL;
    cram_index *e = NULL;

    if (!cidx || !iter)
        return NULL;

    iter->is_cram = 1;
    iter->read_rest = 0;
    iter->off = NULL;
    iter->n_off = 0;
    iter->curr_off = 0;
    iter->i = -1;

    for (i=0; i<iter->n_reg; i++) {

        curr_reg = &iter->reg_list[i];
        tid = curr_reg->tid;

        if (tid >= 0) {
            off = (hts_pair64_max_t*)realloc(off, (n_off + curr_reg->count) * sizeof(hts_pair64_max_t));
            if (!off)
                return NULL;

            for (j=0; j < curr_reg->count; j++) {
                curr_intv = &curr_reg->intervals[j];
                if (curr_intv->end < curr_intv->beg)
                    continue;

                beg = curr_intv->beg;
                end = curr_intv->end;

/* First, fetch the container overlapping 'beg' and assign its file offset to u, then
 * find the container overlapping 'end' and assing the relative end of the slice to v.
 * The cram_ptell function will adjust with the container offset, which is not stored
 * in the index.
 */
                e = cram_index_query(cidx->cram, tid, beg+1, NULL);
                if (e) {
                    off[n_off].u = e->offset;

                    if (end == INT_MAX) {
                       e = cram_index_last(cidx->cram, tid, NULL);
                    } else {
                       e = cram_index_query(cidx->cram, tid, end+1, NULL);
                    }

                    if (e) {
                        off[n_off].v = e->offset + e->slice + e->len;
                        off[n_off].max = (uint64_t)tid<<32 | end;
                        n_off++;
                    } else {
                        hts_log_warning("Could not set offset end for region %d(%s):%d-%d. Skipping", tid, curr_reg->reg, beg, end);
                    }
                } else {
                    hts_log_warning("No index entry for region %d:%d-%d", tid, beg, end);
                }
            }
        } else {
            switch (tid) {
                case HTS_IDX_NOCOOR:
                    e = cram_index_query(cidx->cram, tid, 1, NULL);
                    if (e) {
                        iter->nocoor = 1;
                        iter->nocoor_off = e->offset;
                    } else {
                        hts_log_warning("No index entry for NOCOOR region");
                    }
                    break;
                case HTS_IDX_START:
                    e = cram_index_query(cidx->cram, tid, 1, NULL);
                    if (e) {
                        iter->read_rest = 1;
                        off = (hts_pair64_max_t*)realloc(off, sizeof(hts_pair64_max_t));
                        off[0].u = e->offset;
                        off[0].v = 0;
                        off[0].max = 0;
                        n_off=1;
                    } else {
                        hts_log_warning("No index entries");
                    }
                    break;
                case HTS_IDX_REST:
                    break;
                case HTS_IDX_NONE:
                    iter->finished = 1;
                    break;
                default:
                    hts_log_error("Query with tid=%d not implemented for CRAM files", tid);
            }
        }
    }

    if (n_off) {
        ks_introsort(_off_max, n_off, off);
        // resolve completely contained adjacent blocks
        for (i = 1, l = 0; i < n_off; ++i) {
            if (off[l].v < off[i].v) {
                off[++l] = off[i];
            } else {
                off[l].max = (off[i].max > off[l].max ? off[i].max : off[l].max);
            }
        }
        n_off = l + 1;
        // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
        for (i = 1; i < n_off; ++i)
            if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
        // merge adjacent blocks
        for (i = 1, l = 0; i < n_off; ++i) {
            if (off[l].v>>16 == off[i].u>>16) {
                off[l].v = off[i].v;
                off[l].max = (off[i].max > off[l].max ? off[i].max : off[l].max);
            } else off[++l] = off[i];
        }
        n_off = l + 1;
        iter->n_off = n_off; iter->off = off;
    }

    if(!n_off && !iter->nocoor)
        iter->finished = 1;

    return iter;
}

void hts_itr_destroy(hts_itr_t *iter)
{
    if (iter) { free(iter->off); free(iter->bins.a); free(iter); }
}

void hts_reglist_free(hts_reglist_t *reglist, int count) {

    int i;
    if(reglist) {
        for (i=0;i<count;i++) {
            if (reglist[i].intervals)
                free(reglist[i].intervals);
        }
        free(reglist);
    }
}

void hts_itr_multi_destroy(hts_itr_multi_t *iter) {

    if (iter) {
        if (iter->reg_list && iter->n_reg)
            hts_reglist_free(iter->reg_list, iter->n_reg);

        if (iter->off && iter->n_off)
            free(iter->off);
        free(iter);
    }
}

static inline long long push_digit(long long i, char c)
{
    // ensure subtraction occurs first, avoiding overflow for >= MAX-48 or so
    int digit = c - '0';
    return 10 * i + digit;
}

long long hts_parse_decimal(const char *str, char **strend, int flags)
{
    long long n = 0;
    int decimals = 0, e = 0, lost = 0;
    char sign = '+', esign = '+';
    const char *s;

    while (isspace_c(*str)) str++;
    s = str;

    if (*s == '+' || *s == '-') sign = *s++;
    while (*s)
        if (isdigit_c(*s)) n = push_digit(n, *s++);
        else if (*s == ',' && (flags & HTS_PARSE_THOUSANDS_SEP)) s++;
        else break;

    if (*s == '.') {
        s++;
        while (isdigit_c(*s)) decimals++, n = push_digit(n, *s++);
    }

    if (*s == 'E' || *s == 'e') {
        s++;
        if (*s == '+' || *s == '-') esign = *s++;
        while (isdigit_c(*s)) e = push_digit(e, *s++);
        if (esign == '-') e = -e;
    }

    e -= decimals;
    while (e > 0) n *= 10, e--;
    while (e < 0) lost += n % 10, n /= 10, e++;

    if (lost > 0) {
        hts_log_warning("Discarding fractional part of %.*s", (int)(s - str), str);
    }

    if (strend) {
        *strend = (char *)s;
    } else if (*s) {
        hts_log_warning("Ignoring unknown characters after %.*s[%s]", (int)(s - str), str, s);
    }

    return (sign == '+')? n : -n;
}

const char *hts_parse_reg(const char *s, int *beg, int *end)
{
    char *hyphen;
    const char *colon = strrchr(s, ':');
    if (colon == NULL) {
        *beg = 0; *end = INT_MAX;
        return s + strlen(s);
    }

    *beg = hts_parse_decimal(colon+1, &hyphen, HTS_PARSE_THOUSANDS_SEP) - 1;
    if (*beg < 0) *beg = 0;

    if (*hyphen == '\0') *end = INT_MAX;
    else if (*hyphen == '-') *end = hts_parse_decimal(hyphen+1, NULL, HTS_PARSE_THOUSANDS_SEP);
    else return NULL;

    if (*beg >= *end) return NULL;
    return colon;
}

hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec)
{
    int tid, beg, end;
    const char *q;

    if (strcmp(reg, ".") == 0)
        return itr_query(idx, HTS_IDX_START, 0, 0, readrec);
    else if (strcmp(reg, "*") == 0)
        return itr_query(idx, HTS_IDX_NOCOOR, 0, 0, readrec);

    q = hts_parse_reg(reg, &beg, &end);
    if (q) {
        char tmp_a[1024], *tmp = tmp_a;
        if (q - reg + 1 > 1024)
            if (!(tmp = malloc(q - reg + 1)))
                return NULL;
        strncpy(tmp, reg, q - reg);
        tmp[q - reg] = 0;
        tid = getid(hdr, tmp);
        if (tmp != tmp_a)
            free(tmp);
    }
    else {
        // not parsable as a region, but possibly a sequence named "foo:a"
        tid = getid(hdr, reg);
        beg = 0; end = INT_MAX;
    }

    if (tid < 0) return NULL;
    return itr_query(idx, tid, beg, end, readrec);
}

hts_itr_multi_t *hts_itr_regions(const hts_idx_t *idx, hts_reglist_t *reglist, int count, hts_name2id_f getid, void *hdr, hts_itr_multi_query_func *itr_specific, hts_readrec_func *readrec, hts_seek_func *seek, hts_tell_func *tell) {

    int i;

    if (!reglist)
        return NULL;

    hts_itr_multi_t *itr = (hts_itr_multi_t*)calloc(1, sizeof(hts_itr_multi_t));
    if (itr) {
        itr->n_reg = count;
        itr->readrec = readrec;
        itr->seek = seek;
        itr->tell = tell;
        itr->reg_list = reglist;
        itr->finished = 0;
        itr->nocoor = 0;


        for (i = 0; i < itr->n_reg; i++) {
            if (!strcmp(itr->reg_list[i].reg, ".")) {
                itr->reg_list[i].tid = HTS_IDX_START;
                continue;
            }

            if (!strcmp(itr->reg_list[i].reg, "*")) {
                itr->reg_list[i].tid = HTS_IDX_NOCOOR;
                continue;
            }

            itr->reg_list[i].tid = getid(hdr, reglist[i].reg);
            if (itr->reg_list[i].tid < 0)
                hts_log_warning("Region '%s' specifies an unknown reference name. Continue anyway", reglist[i].reg);
        }

        qsort(itr->reg_list, itr->n_reg, sizeof(hts_reglist_t), compare_regions);
        itr_specific(idx, itr);
    }
    return itr;
}

int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)
{
    int ret, tid, beg, end;
    if (iter == NULL || iter->finished) return -1;
    if (iter->read_rest) {
        if (iter->curr_off) { // seek to the start
            if (bgzf_seek(fp, iter->curr_off, SEEK_SET) < 0) return -1;
            iter->curr_off = 0; // only seek once
        }
        ret = iter->readrec(fp, data, r, &tid, &beg, &end);
        if (ret < 0) iter->finished = 1;
        iter->curr_tid = tid;
        iter->curr_beg = beg;
        iter->curr_end = end;
        return ret;
    }
    // A NULL iter->off should always be accompanied by iter->finished.
    assert(iter->off != NULL);
    for (;;) {
        if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
            if (iter->i == iter->n_off - 1) { ret = -1; break; } // no more chunks
            if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
                if (bgzf_seek(fp, iter->off[iter->i+1].u, SEEK_SET) < 0) return -1;
                iter->curr_off = bgzf_tell(fp);
            }
            ++iter->i;
        }
        if ((ret = iter->readrec(fp, data, r, &tid, &beg, &end)) >= 0) {
            iter->curr_off = bgzf_tell(fp);
            if (tid != iter->tid || beg >= iter->end) { // no need to proceed
                ret = -1; break;
            } else if (end > iter->beg && iter->end > beg) {
                iter->curr_tid = tid;
                iter->curr_beg = beg;
                iter->curr_end = end;
                return ret;
            }
        } else break; // end of file or error
    }
    iter->finished = 1;
    return ret;
}

int hts_itr_multi_next(htsFile *fd, hts_itr_multi_t *iter, void *r)
{
    void *fp;
    int ret, tid, beg, end, i, cr, ci;
    hts_reglist_t *found_reg;

    if (iter == NULL || iter->finished) return -1;

    if (iter->is_cram) {
        fp = fd->fp.cram;
    } else {
        fp = fd->fp.bgzf;
    }

    if (iter->read_rest) {
        if (iter->curr_off) { // seek to the start
            if (iter->seek(fp, iter->curr_off, SEEK_SET) < 0) {
                return -1;
            }
            iter->curr_off = 0; // only seek once
        }

        ret = iter->readrec(fp, fd, r, &tid, &beg, &end);
        if (ret < 0) {
            iter->finished = 1;
        }

        iter->curr_tid = tid;
        iter->curr_beg = beg;
        iter->curr_end = end;

        return ret;
    }
    // A NULL iter->off should always be accompanied by iter->finished.
    assert(iter->off != NULL || iter->nocoor != 0);

    for (;;) {
        if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
            if (iter->i == iter->n_off - 1) { // no more chunks, except NOCOORs
               if (iter->nocoor) {
                   iter->read_rest = 1;
                   iter->curr_off = iter->nocoor_off;

                   return hts_itr_multi_next(fd, iter, r);
               } else {
                   ret = -1; break;
               }
            }

            if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
                if (iter->seek(fp, iter->off[iter->i+1].u, SEEK_SET) < 0) {
                    return -1;
                }

                iter->curr_off = iter->tell(fp);
            }
            ++iter->i;
        }

        ret = iter->readrec(fp, fd, r, &tid, &beg, &end);
        if (ret < 0)
            break;

        iter->curr_off = iter->tell(fp);
        if (tid != iter->curr_tid) {
            hts_reglist_t key;
            key.tid = tid;

            found_reg = (hts_reglist_t *)bsearch(&key, iter->reg_list, iter->n_reg, sizeof(hts_reglist_t), compare_regions);
            if (!found_reg)
                continue;

            iter->curr_reg = (found_reg - iter->reg_list);
            iter->curr_tid = tid;
            iter->curr_intv = 0;
        }

        cr = iter->curr_reg;
        ci = iter->curr_intv;

        if (beg >  iter->off[iter->i].max) {
            iter->curr_off = iter->off[iter->i].v;
            continue;
        }
        if (beg >  iter->reg_list[cr].max_end)
            continue;

        for (i = ci; i < iter->reg_list[cr].count; i++) {
            if (end > iter->reg_list[cr].intervals[i].beg && iter->reg_list[cr].intervals[i].end > beg) {
                iter->curr_beg = beg;
                iter->curr_end = end;
                iter->curr_intv = i;

                return ret;
            }
        }
    }
    iter->finished = 1;

    return ret;
}

/**********************
 *** Retrieve index ***
 **********************/
// Returns -1 if index couldn't be opened.
//         -2 on other errors
static int test_and_fetch(const char *fn, const char **local_fn)
{
    hFILE *remote_hfp;
    FILE *local_fp = NULL;
    uint8_t *buf = NULL;
    int save_errno;

    if (hisremote(fn)) {
        const int buf_size = 1 * 1024 * 1024;
        int l;
        const char *p;
        for (p = fn + strlen(fn) - 1; p >= fn; --p)
            if (*p == '/') break;
        ++p; // p now points to the local file name
        // Attempt to open local file first
        if ((local_fp = fopen((char*)p, "rb")) != 0)
        {
            fclose(local_fp);
            *local_fn = p;
            return 0;
        }
        // Attempt to open remote file. Stay quiet on failure, it is OK to fail when trying first .csi then .tbi index.
        if ((remote_hfp = hopen(fn, "r")) == 0) return -1;
        if ((local_fp = fopen(p, "w")) == 0) {
            hts_log_error("Failed to create file %s in the working directory", p);
            goto fail;
        }
        hts_log_info("Downloading file %s to local directory", fn);
        buf = (uint8_t*)calloc(buf_size, 1);
        if (!buf) {
            hts_log_error("%s", strerror(errno));
            goto fail;
        }
        while ((l = hread(remote_hfp, buf, buf_size)) > 0) {
            if (fwrite(buf, 1, l, local_fp) != l) {
                hts_log_error("Failed to write data to %s : %s",
                              fn, strerror(errno));
                goto fail;
            }
        }
        free(buf);
        if (fclose(local_fp) < 0) {
            hts_log_error("Error closing %s : %s", fn, strerror(errno));
            local_fp = NULL;
            goto fail;
        }
        if (hclose(remote_hfp) != 0) {
            hts_log_error("Failed to close remote file %s", fn);
        }
        *local_fn = p;
        return 0;
    } else {
        hFILE *local_hfp;
        if ((local_hfp = hopen(fn, "r")) == 0) return -1;
        hclose_abruptly(local_hfp);
        *local_fn = fn;
        return 0;
    }

 fail:
    save_errno = errno;
    hclose_abruptly(remote_hfp);
    if (local_fp) fclose(local_fp);
    free(buf);
    errno = save_errno;
    return -2;
}

char *hts_idx_getfn(const char *fn, const char *ext)
{
    int i, l_fn, l_ext, ret;
    char *fnidx;
    const char *local_fn = NULL;
    l_fn = strlen(fn); l_ext = strlen(ext);
    fnidx = (char*)calloc(l_fn + l_ext + 1, 1);
    if (!fnidx) return NULL;
    // First try : append `ext` to `fn`
    strcpy(fnidx, fn); strcpy(fnidx + l_fn, ext);
    if ((ret = test_and_fetch(fnidx, &local_fn)) == -1) {
        // Second try : replace suffix of `fn` with `ext`
        for (i = l_fn - 1; i > 0; --i)
            if (fnidx[i] == '.' || fnidx[i] == '/') break;
        if (fnidx[i] == '.') {
            strcpy(fnidx + i, ext);
            ret = test_and_fetch(fnidx, &local_fn);
        }
    }
    if (ret < 0) {
        free(fnidx);
        return NULL;
    }
    l_fn = strlen(local_fn);
    memmove(fnidx, local_fn, l_fn + 1);
    return fnidx;
}

hts_idx_t *hts_idx_load(const char *fn, int fmt)
{
    char *fnidx;
    hts_idx_t *idx;
    fnidx = hts_idx_getfn(fn, ".csi");
    if (! fnidx) fnidx = hts_idx_getfn(fn, fmt == HTS_FMT_BAI? ".bai" : ".tbi");
    if (fnidx == 0) return 0;

    idx = hts_idx_load2(fn, fnidx);
    free(fnidx);
    return idx;
}

hts_idx_t *hts_idx_load2(const char *fn, const char *fnidx)
{
    // Check that the index file is up to date, the main file might have changed
    struct stat stat_idx,stat_main;
    if ( !stat(fn, &stat_main) && !stat(fnidx, &stat_idx) )
    {
        if ( stat_idx.st_mtime < stat_main.st_mtime )
            hts_log_warning("The index file is older than the data file: %s", fnidx);
    }

    return hts_idx_load_local(fnidx);
}



/**********************
 ***     Memory     ***
 **********************/

/* For use with hts_expand macros *only* */
size_t hts_realloc_or_die(size_t n, size_t m, size_t m_sz, size_t size,
                          int clear, void **ptr, const char *func) {
    /* If new_m and size are both below this limit, multiplying them
       together can't overflow */
    const size_t safe = (size_t) 1 << (sizeof(size_t) * 4);
    void *new_ptr;
    size_t bytes, new_m;

    new_m = n;
    kroundup_size_t(new_m);

    bytes = size * new_m;

    /* Check for overflow.  Both ensure that new_m will fit in m (we make the
       pessimistic assumption that m is signed), and that bytes has not
       wrapped around. */
    if (new_m > (((size_t) 1 << (m_sz * 8 - 1)) - 1)
        || ((size > safe || new_m > safe)
            && bytes / new_m != size)) {
        errno = ENOMEM;
        goto die;
    }

    new_ptr = realloc(*ptr, bytes);
    if (new_ptr == NULL) goto die;

    if (clear) {
        if (new_m > m) {
            memset((char *) new_ptr + m * size, 0, (new_m - m) * size);
        }
    }

    *ptr = new_ptr;

    return new_m;

 die:
    hts_log_error("%s", strerror(errno));
    exit(1);
}

void hts_set_log_level(enum htsLogLevel level)
{
    hts_verbose = level;
}

enum htsLogLevel hts_get_log_level()
{
    return hts_verbose;
}

static char get_severity_tag(enum htsLogLevel severity)
{
    switch (severity) {
    case HTS_LOG_ERROR:
        return 'E';
    case HTS_LOG_WARNING:
        return 'W';
    case HTS_LOG_INFO:
        return 'I';
    case HTS_LOG_DEBUG:
        return 'D';
    case HTS_LOG_TRACE:
        return 'T';
    default:
        break;
    }

    return '*';
}

void hts_log(enum htsLogLevel severity, const char *context, const char *format, ...)
{
    int save_errno = errno;
    if (severity <= hts_verbose) {
        va_list argptr;

        fprintf(stderr, "[%c::%s] ", get_severity_tag(severity), context);

        va_start(argptr, format);
        vfprintf(stderr, format, argptr);
        va_end(argptr);

        fprintf(stderr, "\n");
    }
    errno = save_errno;
}
