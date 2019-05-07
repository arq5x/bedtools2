/// @file htslib/hts.h
/// Format-neutral I/O, indexing, and iterator API functions.
/*
    Copyright (C) 2012-2016 Genome Research Ltd.
    Copyright (C) 2010, 2012 Broad Institute.
    Portions copyright (C) 2003-2006, 2008-2010 by Heng Li <lh3@live.co.uk>

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

#ifndef HTSLIB_HTS_H
#define HTSLIB_HTS_H

#include <stddef.h>
#include <stdint.h>

#include "hts_defs.h"
#include "hts_log.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HTS_BGZF_TYPEDEF
typedef struct BGZF BGZF;
#define HTS_BGZF_TYPEDEF
#endif
struct cram_fd;
struct hFILE;
struct hFILE_callback_ops;
struct hts_tpool;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/**
 * @hideinitializer
 * Macro to expand a dynamic array of a given type
 *
 * @param         type_t The type of the array elements
 * @param[in]     n      Requested number of elements of type type_t
 * @param[in,out] m      Size of memory allocated
 * @param[in,out] ptr    Pointer to the array
 *
 * @discussion
 * The array *ptr will be expanded if necessary so that it can hold @p n
 * or more elements.  If the array is expanded then the new size will be
 * written to @p m and the value in @ptr may change.
 *
 * It must be possible to take the address of @p ptr and @p m must be usable
 * as an lvalue.
 *
 * @bug
 * If the memory allocation fails, this will call exit(1).  This is
 * not ideal behaviour in a library.
 */
#define hts_expand(type_t, n, m, ptr) do {                              \
        if ((n) > (m)) {                                                \
            size_t hts_realloc_or_die(size_t, size_t, size_t, size_t,   \
                                      int, void **, const char *);      \
            (m) = hts_realloc_or_die((n) >= 1 ? (n) : 1, (m), sizeof(m), \
                                     sizeof(type_t),  0,                \
                                     (void **)&(ptr), __func__);        \
        }                                                               \
    } while (0)

/**
 * @hideinitializer
 * Macro to expand a dynamic array, zeroing any newly-allocated memory
 *
 * @param         type_t The type of the array elements
 * @param[in]     n      Requested number of elements of type type_t
 * @param[in,out] m      Size of memory allocated
 * @param[in,out] ptr    Pointer to the array
 *
 * @discussion
 * As for hts_expand(), except the bytes that make up the array elements
 * between the old and new values of @p m are set to zero using memset().
 *
 * @bug
 * If the memory allocation fails, this will call exit(1).  This is
 * not ideal behaviour in a library.
 */


#define hts_expand0(type_t, n, m, ptr) do {                             \
        if ((n) > (m)) {                                                \
            size_t hts_realloc_or_die(size_t, size_t, size_t, size_t,   \
                                      int, void **, const char *);      \
            (m) = hts_realloc_or_die((n) >= 1 ? (n) : 1, (m), sizeof(m), \
                                     sizeof(type_t), 1,                 \
                                     (void **)&(ptr), __func__);        \
        }                                                               \
    } while (0)

/************
 * File I/O *
 ************/

// Add new entries only at the end (but before the *_maximum entry)
// of these enums, as their numbering is part of the htslib ABI.

enum htsFormatCategory {
    unknown_category,
    sequence_data,    // Sequence data -- SAM, BAM, CRAM, etc
    variant_data,     // Variant calling data -- VCF, BCF, etc
    index_file,       // Index file associated with some data file
    region_list,      // Coordinate intervals or regions -- BED, etc
    category_maximum = 32767
};

enum htsExactFormat {
    unknown_format,
    binary_format, text_format,
    sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed,
    htsget,
    json HTS_DEPRECATED_ENUM("Use htsExactFormat 'htsget' instead") = htsget,
    format_maximum = 32767
};

enum htsCompression {
    no_compression, gzip, bgzf, custom,
    compression_maximum = 32767
};

typedef struct htsFormat {
    enum htsFormatCategory category;
    enum htsExactFormat format;
    struct { short major, minor; } version;
    enum htsCompression compression;
    short compression_level;  // currently unused
    void *specific;  // format specific options; see struct hts_opt.
} htsFormat;

// Maintainers note htsFile cannot be an opaque structure because some of its
// fields are part of libhts.so's ABI (hence these fields must not be moved):
//  - fp is used in the public sam_itr_next()/etc macros
//  - is_bin is used directly in samtools <= 1.1 and bcftools <= 1.1
//  - is_write and is_cram are used directly in samtools <= 1.1
//  - fp is used directly in samtools (up to and including current develop)
//  - line is used directly in bcftools (up to and including current develop)
typedef struct {
    uint32_t is_bin:1, is_write:1, is_be:1, is_cram:1, is_bgzf:1, dummy:27;
    int64_t lineno;
    kstring_t line;
    char *fn, *fn_aux;
    union {
        BGZF *bgzf;
        struct cram_fd *cram;
        struct hFILE *hfile;
    } fp;
    htsFormat format;
} htsFile;

// A combined thread pool and queue allocation size.
// The pool should already be defined, but qsize may be zero to
// indicate an appropriate queue size is taken from the pool.
//
// Reasons for explicitly setting it could be where many more file
// descriptors are in use than threads, so keeping memory low is
// important.
typedef struct {
    struct hts_tpool *pool; // The shared thread pool itself
    int qsize;    // Size of I/O queue to use for this fp
} htsThreadPool;

// REQUIRED_FIELDS
enum sam_fields {
    SAM_QNAME = 0x00000001,
    SAM_FLAG  = 0x00000002,
    SAM_RNAME = 0x00000004,
    SAM_POS   = 0x00000008,
    SAM_MAPQ  = 0x00000010,
    SAM_CIGAR = 0x00000020,
    SAM_RNEXT = 0x00000040,
    SAM_PNEXT = 0x00000080,
    SAM_TLEN  = 0x00000100,
    SAM_SEQ   = 0x00000200,
    SAM_QUAL  = 0x00000400,
    SAM_AUX   = 0x00000800,
    SAM_RGAUX = 0x00001000,
};

// Mostly CRAM only, but this could also include other format options
enum hts_fmt_option {
    // CRAM specific
    CRAM_OPT_DECODE_MD,
    CRAM_OPT_PREFIX,
    CRAM_OPT_VERBOSITY,  // obsolete, use hts_set_log_level() instead
    CRAM_OPT_SEQS_PER_SLICE,
    CRAM_OPT_SLICES_PER_CONTAINER,
    CRAM_OPT_RANGE,
    CRAM_OPT_VERSION,    // rename to cram_version?
    CRAM_OPT_EMBED_REF,
    CRAM_OPT_IGNORE_MD5,
    CRAM_OPT_REFERENCE,  // make general
    CRAM_OPT_MULTI_SEQ_PER_SLICE,
    CRAM_OPT_NO_REF,
    CRAM_OPT_USE_BZIP2,
    CRAM_OPT_SHARED_REF,
    CRAM_OPT_NTHREADS,   // deprecated, use HTS_OPT_NTHREADS
    CRAM_OPT_THREAD_POOL,// make general
    CRAM_OPT_USE_LZMA,
    CRAM_OPT_USE_RANS,
    CRAM_OPT_REQUIRED_FIELDS,
    CRAM_OPT_LOSSY_NAMES,
    CRAM_OPT_BASES_PER_SLICE,
    CRAM_OPT_STORE_MD,
    CRAM_OPT_STORE_NM,

    // General purpose
    HTS_OPT_COMPRESSION_LEVEL = 100,
    HTS_OPT_NTHREADS,
    HTS_OPT_THREAD_POOL,
    HTS_OPT_CACHE_SIZE,
    HTS_OPT_BLOCK_SIZE,
};

// For backwards compatibility
#define cram_option hts_fmt_option

typedef struct hts_opt {
    char *arg;                // string form, strdup()ed
    enum hts_fmt_option opt;  // tokenised key
    union {                   // ... and value
        int i;
        char *s;
    } val;
    struct hts_opt *next;
} hts_opt;

#define HTS_FILE_OPTS_INIT {{0},0}

/**********************
 * Exported functions *
 **********************/

/*
 * Parses arg and appends it to the option list.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int hts_opt_add(hts_opt **opts, const char *c_arg);

/*
 * Applies an hts_opt option list to a given htsFile.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int hts_opt_apply(htsFile *fp, hts_opt *opts);

/*
 * Frees an hts_opt list.
 */
void hts_opt_free(hts_opt *opts);

/*
 * Accepts a string file format (sam, bam, cram, vcf, bam) optionally
 * followed by a comma separated list of key=value options and splits
 * these up into the fields of htsFormat struct.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
int hts_parse_format(htsFormat *opt, const char *str);

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
int hts_parse_opt_list(htsFormat *opt, const char *str);

/*! @abstract Table for converting a nucleotide character to 4-bit encoding.
The input character may be either an IUPAC ambiguity code, '=' for 0, or
'0'/'1'/'2'/'3' for a result of 1/2/4/8.  The result is encoded as 1/2/4/8
for A/C/G/T or combinations of these bits for ambiguous bases.
*/
extern const unsigned char seq_nt16_table[256];

/*! @abstract Table for converting a 4-bit encoded nucleotide to an IUPAC
ambiguity code letter (or '=' when given 0).
*/
extern const char seq_nt16_str[];

/*! @abstract Table for converting a 4-bit encoded nucleotide to about 2 bits.
Returns 0/1/2/3 for 1/2/4/8 (i.e., A/C/G/T), or 4 otherwise (0 or ambiguous).
*/
extern const int seq_nt16_int[];

/*!
  @abstract  Get the htslib version number
  @return    For released versions, a string like "N.N[.N]"; or git describe
  output if using a library built within a Git repository.
*/
const char *hts_version(void);

/*!
  @abstract    Determine format by peeking at the start of a file
  @param fp    File opened for reading, positioned at the beginning
  @param fmt   Format structure that will be filled out on return
  @return      0 for success, or negative if an error occurred.
*/
int hts_detect_format(struct hFILE *fp, htsFormat *fmt);

/*!
  @abstract    Get a human-readable description of the file format
  @param fmt   Format structure holding type, version, compression, etc.
  @return      Description string, to be freed by the caller after use.
*/
char *hts_format_description(const htsFormat *format);

/*!
  @abstract       Open a SAM/BAM/CRAM/VCF/BCF/etc file
  @param fn       The file name or "-" for stdin/stdout
  @param mode     Mode matching / [rwa][bceguxz0-9]* /
  @discussion
      With 'r' opens for reading; any further format mode letters are ignored
      as the format is detected by checking the first few bytes or BGZF blocks
      of the file.  With 'w' or 'a' opens for writing or appending, with format
      specifier letters:
        b  binary format (BAM, BCF, etc) rather than text (SAM, VCF, etc)
        c  CRAM format
        g  gzip compressed
        u  uncompressed
        z  bgzf compressed
        [0-9]  zlib compression level
      and with non-format option letters (for any of 'r'/'w'/'a'):
        e  close the file on exec(2) (opens with O_CLOEXEC, where supported)
        x  create the file exclusively (opens with O_EXCL, where supported)
      Note that there is a distinction between 'u' and '0': the first yields
      plain uncompressed output whereas the latter outputs uncompressed data
      wrapped in the zlib format.
  @example
      [rw]b  .. compressed BCF, BAM, FAI
      [rw]bu .. uncompressed BCF
      [rw]z  .. compressed VCF
      [rw]   .. uncompressed VCF
*/
htsFile *hts_open(const char *fn, const char *mode);

/*!
  @abstract      Open a SAM/BAM/CRAM/VCF/BCF/etc from stream callbacks
  @param ops     The IO operation callback
  @param mode    The mode string
  @discussion    
     See hts_open for description of mode string.
	 This function is useful when data sources other than file/pipe
	 should be consumed.
*/
htsFile *hts_open_callback(const char* fn, struct hFILE_callback_ops* ops, const char* mode);

/*!
  @abstract       Open a SAM/BAM/CRAM/VCF/BCF/etc file
  @param fn       The file name or "-" for stdin/stdout
  @param mode     Open mode, as per hts_open()
  @param fmt      Optional format specific parameters
  @discussion
      See hts_open() for description of fn and mode.
      // TODO Update documentation for s/opts/fmt/
      Opts contains a format string (sam, bam, cram, vcf, bcf) which will,
      if defined, override mode.  Opts also contains a linked list of hts_opt
      structures to apply to the open file handle.  These can contain things
      like pointers to the reference or information on compression levels,
      block sizes, etc.
*/
htsFile *hts_open_format(const char *fn, const char *mode, const htsFormat *fmt);

/*!
  @abstract       Open an existing stream as a SAM/BAM/CRAM/VCF/BCF/etc file
  @param fn       The already-open file handle
  @param mode     Open mode, as per hts_open()
*/
htsFile *hts_hopen(struct hFILE *fp, const char *fn, const char *mode);

/*!
  @abstract  Close a file handle, flushing buffered data for output streams
  @param fp  The file handle to be closed
  @return    0 for success, or negative if an error occurred.
*/
int hts_close(htsFile *fp);

/*!
  @abstract  Returns the file's format information
  @param fp  The file handle
  @return    Read-only pointer to the file's htsFormat.
*/
const htsFormat *hts_get_format(htsFile *fp);

/*!
  @ abstract      Returns a string containing the file format extension.
  @ param format  Format structure containing the file type.
  @ return        A string ("sam", "bam", etc) or "?" for unknown formats.
 */
const char *hts_format_file_extension(const htsFormat *format);

/*!
  @abstract  Sets a specified CRAM option on the open file handle.
  @param fp  The file handle open the open file.
  @param opt The CRAM_OPT_* option.
  @param ... Optional arguments, dependent on the option used.
  @return    0 for success, or negative if an error occurred.
*/
int hts_set_opt(htsFile *fp, enum hts_fmt_option opt, ...);

int hts_getline(htsFile *fp, int delimiter, kstring_t *str);
char **hts_readlines(const char *fn, int *_n);
/*!
    @abstract       Parse comma-separated list or read list from a file
    @param list     File name or comma-separated list
    @param is_file
    @param _n       Size of the output array (number of items read)
    @return         NULL on failure or pointer to newly allocated array of
                    strings
*/
char **hts_readlist(const char *fn, int is_file, int *_n);

/*!
  @abstract  Create extra threads to aid compress/decompression for this file
  @param fp  The file handle
  @param n   The number of worker threads to create
  @return    0 for success, or negative if an error occurred.
  @notes     This function creates non-shared threads for use solely by fp.
             The hts_set_thread_pool function is the recommended alternative.
*/
int hts_set_threads(htsFile *fp, int n);

/*!
  @abstract  Create extra threads to aid compress/decompression for this file
  @param fp  The file handle
  @param p   A pool of worker threads, previously allocated by hts_create_threads().
  @return    0 for success, or negative if an error occurred.
*/
int hts_set_thread_pool(htsFile *fp, htsThreadPool *p);

/*!
  @abstract  Adds a cache of decompressed blocks, potentially speeding up seeks.
             This may not work for all file types (currently it is bgzf only).
  @param fp  The file handle
  @param n   The size of cache, in bytes
*/
void hts_set_cache_size(htsFile *fp, int n);

/*!
  @abstract  Set .fai filename for a file opened for reading
  @return    0 for success, negative on failure
  @discussion
      Called before *_hdr_read(), this provides the name of a .fai file
      used to provide a reference list if the htsFile contains no @SQ headers.
*/
int hts_set_fai_filename(htsFile *fp, const char *fn_aux);


/*!
  @abstract  Determine whether a given htsFile contains a valid EOF block
  @return    3 for a non-EOF checkable filetype;
             2 for an unseekable file type where EOF cannot be checked;
             1 for a valid EOF block;
             0 for if the EOF marker is absent when it should be present;
            -1 (with errno set) on failure
  @discussion
      Check if the BGZF end-of-file (EOF) marker is present
*/
int hts_check_EOF(htsFile *fp);

/************
 * Indexing *
 ************/

/*!
These HTS_IDX_* macros are used as special tid values for hts_itr_query()/etc,
producing iterators operating as follows:
 - HTS_IDX_NOCOOR iterates over unmapped reads sorted at the end of the file
 - HTS_IDX_START  iterates over the entire file
 - HTS_IDX_REST   iterates from the current position to the end of the file
 - HTS_IDX_NONE   always returns "no more alignment records"
When one of these special tid values is used, beg and end are ignored.
When REST or NONE is used, idx is also ignored and may be NULL.
*/
#define HTS_IDX_NOCOOR (-2)
#define HTS_IDX_START  (-3)
#define HTS_IDX_REST   (-4)
#define HTS_IDX_NONE   (-5)

#define HTS_FMT_CSI 0
#define HTS_FMT_BAI 1
#define HTS_FMT_TBI 2
#define HTS_FMT_CRAI 3

struct __hts_idx_t;
typedef struct __hts_idx_t hts_idx_t;

typedef struct {
    uint32_t beg, end;
} hts_pair32_t;

typedef struct {
    uint64_t u, v;
} hts_pair64_t;

typedef struct {
    uint64_t u, v;
    uint64_t max;
} hts_pair64_max_t;

typedef struct {
    const char *reg;
    int tid;
    hts_pair32_t *intervals;
    uint32_t count;
    uint32_t min_beg, max_end;
} hts_reglist_t;

typedef int hts_readrec_func(BGZF *fp, void *data, void *r, int *tid, int *beg, int *end);
typedef int hts_seek_func(void *fp, int64_t offset, int where);
typedef int64_t hts_tell_func(void *fp);

typedef struct {
    uint32_t read_rest:1, finished:1, is_cram:1, dummy:29;
    int tid, beg, end, n_off, i;
    int curr_tid, curr_beg, curr_end;
    uint64_t curr_off;
    hts_pair64_t *off;
    hts_readrec_func *readrec;
    struct {
        int n, m;
        int *a;
    } bins;
} hts_itr_t;

typedef struct {
    int key;
    uint64_t min_off, max_off;
} aux_key_t;

typedef struct {
    uint32_t read_rest:1, finished:1, is_cram:1, nocoor:1, dummy:28;
    hts_reglist_t *reg_list;
    int n_reg, i;
    int curr_tid, curr_intv, curr_beg, curr_end, curr_reg;
    hts_pair64_max_t *off;
    int n_off;
    uint64_t curr_off, nocoor_off;
    hts_readrec_func *readrec;
    hts_seek_func *seek;
    hts_tell_func *tell;
} hts_itr_multi_t;

    #define hts_bin_first(l) (((1<<(((l)<<1) + (l))) - 1) / 7)
    #define hts_bin_parent(l) (((l) - 1) >> 3)

    hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls);
    void hts_idx_destroy(hts_idx_t *idx);
    int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped);
    void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset);

/// Save an index to a file
/** @param idx  Index to be written
    @param fn   Input BAM/BCF/etc filename, to which .bai/.csi/etc will be added
    @param fmt  One of the HTS_FMT_* index formats
    @return  0 if successful, or negative if an error occurred.
*/
int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt) HTS_RESULT_USED;

/// Save an index to a specific file
/** @param idx    Index to be written
    @param fn     Input BAM/BCF/etc filename
    @param fnidx  Output filename, or NULL to add .bai/.csi/etc to @a fn
    @param fmt    One of the HTS_FMT_* index formats
    @return  0 if successful, or negative if an error occurred.
*/
int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt) HTS_RESULT_USED;

/// Load an index file
/** @param fn   BAM/BCF/etc filename, to which .bai/.csi/etc will be added or
                the extension substituted, to search for an existing index file
    @param fmt  One of the HTS_FMT_* index formats
    @return  The index, or NULL if an error occurred.
*/
hts_idx_t *hts_idx_load(const char *fn, int fmt);

/// Load a specific index file
/** @param fn     Input BAM/BCF/etc filename
    @param fnidx  The input index filename
    @return  The index, or NULL if an error occurred.
*/
hts_idx_t *hts_idx_load2(const char *fn, const char *fnidx);


/// Get extra index meta-data
/** @param idx    The index
    @param l_meta Pointer to where the length of the extra data is stored
    @return Pointer to the extra data if present; NULL otherwise

    Indexes (both .tbi and .csi) made by tabix include extra data about
    the indexed file.  The returns a pointer to this data.  Note that the
    data is stored exactly as it is in the index.  Callers need to interpret
    the results themselves, including knowing what sort of data to expect;
    byte swapping etc.
*/
uint8_t *hts_idx_get_meta(hts_idx_t *idx, uint32_t *l_meta);

/// Set extra index meta-data
/** @param idx     The index
    @param l_meta  Length of data
    @param meta    Pointer to the extra data
    @param is_copy If not zero, a copy of the data is taken
    @return 0 on success; -1 on failure (out of memory).

    Sets the data that is returned by hts_idx_get_meta().

    If is_copy != 0, a copy of the input data is taken.  If not, ownership of
    the data pointed to by *meta passes to the index.
*/
int hts_idx_set_meta(hts_idx_t *idx, uint32_t l_meta, uint8_t *meta, int is_copy);

    int hts_idx_get_stat(const hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped);
    uint64_t hts_idx_get_n_no_coor(const hts_idx_t* idx);


#define HTS_PARSE_THOUSANDS_SEP 1  ///< Ignore ',' separators within numbers

/// Parse a numeric string
/** The number may be expressed in scientific notation, and optionally may
    contain commas in the integer part (before any decimal point or E notation).
    @param str     String to be parsed
    @param strend  If non-NULL, set on return to point to the first character
                   in @a str after those forming the parsed number
    @param flags   Or'ed-together combination of HTS_PARSE_* flags
    @return  Converted value of the parsed number.

    When @a strend is NULL, a warning will be printed (if hts_verbose is HTS_LOG_WARNING
    or more) if there are any trailing characters after the number.
*/
long long hts_parse_decimal(const char *str, char **strend, int flags);

/// Parse a "CHR:START-END"-style region string
/** @param str  String to be parsed
    @param beg  Set on return to the 0-based start of the region
    @param end  Set on return to the 1-based end of the region
    @return  Pointer to the colon or '\0' after the reference sequence name,
             or NULL if @a str could not be parsed.
*/
const char *hts_parse_reg(const char *str, int *beg, int *end);

    hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec);
    void hts_itr_destroy(hts_itr_t *iter);

    typedef int (*hts_name2id_f)(void*, const char*);
    typedef const char *(*hts_id2name_f)(void*, int);
    typedef hts_itr_t *hts_itr_query_func(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec);

    hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec);
    int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data) HTS_RESULT_USED;
    const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr); // free only the array, not the values

/**********************************
 * Iterator with multiple regions *
 **********************************/

typedef hts_itr_multi_t *hts_itr_multi_query_func(const hts_idx_t *idx, hts_itr_multi_t *itr);
hts_itr_multi_t *hts_itr_multi_bam(const hts_idx_t *idx, hts_itr_multi_t *iter);
hts_itr_multi_t *hts_itr_multi_cram(const hts_idx_t *idx, hts_itr_multi_t *iter);
hts_itr_multi_t *hts_itr_regions(const hts_idx_t *idx, hts_reglist_t *reglist, int count, hts_name2id_f getid, void *hdr, hts_itr_multi_query_func *itr_specific, hts_readrec_func *readrec, hts_seek_func *seek, hts_tell_func *tell);
int hts_itr_multi_next(htsFile *fd, hts_itr_multi_t *iter, void *r);
void hts_reglist_free(hts_reglist_t *reglist, int count);
void hts_itr_multi_destroy(hts_itr_multi_t *iter);


    /**
     * hts_file_type() - Convenience function to determine file type
     * DEPRECATED:  This function has been replaced by hts_detect_format().
     * It and these FT_* macros will be removed in a future HTSlib release.
     */
    #define FT_UNKN   0
    #define FT_GZ     1
    #define FT_VCF    2
    #define FT_VCF_GZ (FT_GZ|FT_VCF)
    #define FT_BCF    (1<<2)
    #define FT_BCF_GZ (FT_GZ|FT_BCF)
    #define FT_STDIN  (1<<3)
    int hts_file_type(const char *fname);


/***************************
 * Revised MAQ error model *
 ***************************/

struct errmod_t;
typedef struct errmod_t errmod_t;

errmod_t *errmod_init(double depcorr);
void errmod_destroy(errmod_t *em);

/*
    n: number of bases
    m: maximum base
    bases[i]: qual:6, strand:1, base:4
    q[i*m+j]: phred-scaled likelihood of (i,j)
 */
int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *q);


/*****************************************************
 * Probabilistic banded glocal alignment             *
 * See https://doi.org/10.1093/bioinformatics/btr076 *
 *****************************************************/

typedef struct probaln_par_t {
    float d, e;
    int bw;
} probaln_par_t;

/// Perform probabilistic banded glocal alignment
/** @param      ref     Reference sequence
    @param      l_ref   Length of reference
    @param      query   Query sequence
    @param      l_query Length of query sequence
    @param      iqual   Query base qualities
    @param      c       Alignment parameters
    @param[out] state   Output alignment
    @param[out] q    Phred scaled posterior probability of state[i] being wrong
    @return     Phred-scaled likelihood score, or INT_MIN on failure.

The reference and query sequences are coded using integers 0,1,2,3,4 for
bases A,C,G,T,N respectively (N here is for any ambiguity code).

On output, state and q are arrays of length l_query. The higher 30
bits give the reference position the query base is matched to and the
lower two bits can be 0 (an alignment match) or 1 (an
insertion). q[i] gives the phred scaled posterior probability of
state[i] being wrong.

On failure, errno will be set to EINVAL if the values of l_ref or l_query
were invalid; or ENOMEM if a memory allocation failed.
*/

int probaln_glocal(const uint8_t *ref, int l_ref, const uint8_t *query, int l_query, const uint8_t *iqual, const probaln_par_t *c, int *state, uint8_t *q);


    /**********************
     * MD5 implementation *
     **********************/

    struct hts_md5_context;
    typedef struct hts_md5_context hts_md5_context;

    /*! @abstract   Intialises an MD5 context.
     *  @discussion
     *    The expected use is to allocate an hts_md5_context using
     *    hts_md5_init().  This pointer is then passed into one or more calls
     *    of hts_md5_update() to compute successive internal portions of the
     *    MD5 sum, which can then be externalised as a full 16-byte MD5sum
     *    calculation by calling hts_md5_final().  This can then be turned
     *    into ASCII via hts_md5_hex().
     *
     *    To dealloate any resources created by hts_md5_init() call the
     *    hts_md5_destroy() function.
     *
     *  @return     hts_md5_context pointer on success, NULL otherwise.
     */
    hts_md5_context *hts_md5_init(void);

    /*! @abstract Updates the context with the MD5 of the data. */
    void hts_md5_update(hts_md5_context *ctx, const void *data, unsigned long size);

    /*! @abstract Computes the final 128-bit MD5 hash from the given context */
    void hts_md5_final(unsigned char *digest, hts_md5_context *ctx);

    /*! @abstract Resets an md5_context to the initial state, as returned
     *            by hts_md5_init().
     */
    void hts_md5_reset(hts_md5_context *ctx);

    /*! @abstract Converts a 128-bit MD5 hash into a 33-byte nul-termninated
     *            hex string.
     */
    void hts_md5_hex(char *hex, const unsigned char *digest);

    /*! @abstract Deallocates any memory allocated by hts_md5_init. */
    void hts_md5_destroy(hts_md5_context *ctx);


static inline int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
{
    int l, s = min_shift, t = ((1<<((n_lvls<<1) + n_lvls)) - 1) / 7;
    for (--end, l = n_lvls; l > 0; --l, s += 3, t -= 1<<((l<<1)+l))
        if (beg>>s == end>>s) return t + (beg>>s);
    return 0;
}

static inline int hts_bin_bot(int bin, int n_lvls)
{
    int l, b;
    for (l = 0, b = bin; b; ++l, b = hts_bin_parent(b)); // compute the level of bin
    return (bin - hts_bin_first(l)) << (n_lvls - l) * 3;
}

/**************
 * Endianness *
 **************/

static inline int ed_is_big(void)
{
    long one= 1;
    return !(*((char *)(&one)));
}
static inline uint16_t ed_swap_2(uint16_t v)
{
    return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *ed_swap_2p(void *x)
{
    *(uint16_t*)x = ed_swap_2(*(uint16_t*)x);
    return x;
}
static inline uint32_t ed_swap_4(uint32_t v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *ed_swap_4p(void *x)
{
    *(uint32_t*)x = ed_swap_4(*(uint32_t*)x);
    return x;
}
static inline uint64_t ed_swap_8(uint64_t v)
{
    v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
    v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
    return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *ed_swap_8p(void *x)
{
    *(uint64_t*)x = ed_swap_8(*(uint64_t*)x);
    return x;
}

#ifdef __cplusplus
}
#endif

#endif
