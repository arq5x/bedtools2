/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013-2017 Genome Research Ltd

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include <sys/types.h>
#include <inttypes.h>

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/thread_pool.h"
#include "htslib/hts_endian.h"
#include "cram/pooled_alloc.h"

#define BGZF_CACHE
#define BGZF_MT

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


/* BGZF/GZIP header (speciallized from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  B   C

  BGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^16 bytes and adds and an extra "BC" field in the gzip header which
  records the size.

*/
static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";

#ifdef BGZF_CACHE
typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
#include "htslib/khash.h"
KHASH_MAP_INIT_INT64(cache, cache_t)
#endif

#ifdef BGZF_MT

typedef struct bgzf_job {
    BGZF *fp;
    unsigned char comp_data[BGZF_MAX_BLOCK_SIZE];
    size_t comp_len;
    unsigned char uncomp_data[BGZF_MAX_BLOCK_SIZE];
    size_t uncomp_len;
    int errcode;
    int64_t block_address;
    int hit_eof;
} bgzf_job;

enum mtaux_cmd {
    NONE = 0,
    SEEK,
    HAS_EOF,
    CLOSE,
};

typedef struct bgzf_mtaux_t {
    // Memory pool for bgzf_job structs, to avoid many malloc/free
    pool_alloc_t *job_pool;
    bgzf_job *curr_job;

    // Thread pool
    int n_threads;
    int own_pool;
    hts_tpool *pool;

    // Output queue holding completed bgzf_jobs
    hts_tpool_process *out_queue;

    // I/O thread.
    pthread_t io_task;
    pthread_mutex_t job_pool_m;
    int jobs_pending; // number of jobs waiting
    int flush_pending;
    void *free_block;
    int hit_eof;  // r/w entirely within main thread

    // Message passing to the reader thread; eg seek requests
    int errcode;
    uint64_t block_address;
    int eof;
    pthread_mutex_t command_m; // Set whenever fp is being updated
    pthread_cond_t command_c;
    enum mtaux_cmd command;
} mtaux_t;
#endif

typedef struct
{
    uint64_t uaddr;  // offset w.r.t. uncompressed data
    uint64_t caddr;  // offset w.r.t. compressed data
}
bgzidx1_t;

struct __bgzidx_t
{
    int noffs, moffs;       // the size of the index, n:used, m:allocated
    bgzidx1_t *offs;        // offsets
    uint64_t ublock_addr;   // offset of the current block (uncompressed data)
};

void bgzf_index_destroy(BGZF *fp);
int bgzf_index_add_block(BGZF *fp);
static void mt_destroy(mtaux_t *mt);

static inline void packInt16(uint8_t *buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}

static inline void packInt32(uint8_t *buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

static const char *bgzf_zerr(int errnum, z_stream *zs)
{
    static char buffer[32];

    /* Return zs->msg if available.
       zlib doesn't set this very reliably.  Looking at the source suggests
       that it may get set to a useful message for deflateInit2, inflateInit2
       and inflate when it returns Z_DATA_ERROR. For inflate with other
       return codes, deflate, deflateEnd and inflateEnd it doesn't appear
       to be useful.  For the likely non-useful cases, the caller should
       pass NULL into zs. */

    if (zs && zs->msg) return zs->msg;

    // gzerror OF((gzFile file, int *errnum)
    switch (errnum) {
    case Z_ERRNO:
        return strerror(errno);
    case Z_STREAM_ERROR:
        return "invalid parameter/compression level, or inconsistent stream state";
    case Z_DATA_ERROR:
        return "invalid or incomplete IO";
    case Z_MEM_ERROR:
        return "out of memory";
    case Z_BUF_ERROR:
        return "progress temporarily not possible, or in() / out() returned an error";
    case Z_VERSION_ERROR:
        return "zlib version mismatch";
    case Z_OK: // 0: maybe gzgets error Z_NULL
    default:
        snprintf(buffer, sizeof(buffer), "[%d] unknown", errnum);
        return buffer;  // FIXME: Not thread-safe.
    }
}

static BGZF *bgzf_read_init(hFILE *hfpr)
{
    BGZF *fp;
    uint8_t magic[18];
    ssize_t n = hpeek(hfpr, magic, 18);
    if (n < 0) return NULL;

    fp = (BGZF*)calloc(1, sizeof(BGZF));
    if (fp == NULL) return NULL;

    fp->is_write = 0;
    fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
    if (fp->uncompressed_block == NULL) { free(fp); return NULL; }
    fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;
    fp->is_compressed = (n==18 && magic[0]==0x1f && magic[1]==0x8b);
    fp->is_gzip = ( !fp->is_compressed || ((magic[3]&4) && memcmp(&magic[12], "BC\2\0",4)==0) ) ? 0 : 1;
#ifdef BGZF_CACHE
    fp->cache = kh_init(cache);
#endif
    return fp;
}

// get the compress level from the mode string: compress_level==-1 for the default level, -2 plain uncompressed
static int mode2level(const char *mode)
{
    int i, compress_level = -1;
    for (i = 0; mode[i]; ++i)
        if (mode[i] >= '0' && mode[i] <= '9') break;
    if (mode[i]) compress_level = (int)mode[i] - '0';
    if (strchr(mode, 'u')) compress_level = -2;
    return compress_level;
}
static BGZF *bgzf_write_init(const char *mode)
{
    BGZF *fp;
    fp = (BGZF*)calloc(1, sizeof(BGZF));
    if (fp == NULL) goto mem_fail;
    fp->is_write = 1;
    int compress_level = mode2level(mode);
    if ( compress_level==-2 )
    {
        fp->is_compressed = 0;
        return fp;
    }
    fp->is_compressed = 1;

    fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
    if (fp->uncompressed_block == NULL) goto mem_fail;
    fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;

    fp->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1
    if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;
    if ( strchr(mode,'g') )
    {
        // gzip output
        fp->is_gzip = 1;
        fp->gz_stream = (z_stream*)calloc(1,sizeof(z_stream));
        if (fp->gz_stream == NULL) goto mem_fail;
        fp->gz_stream->zalloc = NULL;
        fp->gz_stream->zfree  = NULL;
        fp->gz_stream->msg    = NULL;

        int ret = deflateInit2(fp->gz_stream, fp->compress_level, Z_DEFLATED, 15|16, 8, Z_DEFAULT_STRATEGY);
        if (ret!=Z_OK) {
            hts_log_error("Call to deflateInit2 failed: %s", bgzf_zerr(ret, fp->gz_stream));
            goto fail;
        }
    }
    return fp;

mem_fail:
    hts_log_error("%s", strerror(errno));

fail:
    if (fp != NULL) {
        free(fp->uncompressed_block);
        free(fp->gz_stream);
        free(fp);
    }
    return NULL;
}

BGZF *bgzf_open(const char *path, const char *mode)
{
    BGZF *fp = 0;
    assert(compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);
    if (strchr(mode, 'r')) {
        hFILE *fpr;
        if ((fpr = hopen(path, mode)) == 0) return 0;
        fp = bgzf_read_init(fpr);
        if (fp == 0) { hclose_abruptly(fpr); return NULL; }
        fp->fp = fpr;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        hFILE *fpw;
        if ((fpw = hopen(path, mode)) == 0) return 0;
        fp = bgzf_write_init(mode);
        if (fp == NULL) return NULL;
        fp->fp = fpw;
    }
    else { errno = EINVAL; return 0; }

    fp->is_be = ed_is_big();
    return fp;
}

BGZF *bgzf_dopen(int fd, const char *mode)
{
    BGZF *fp = 0;
    assert(compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);
    if (strchr(mode, 'r')) {
        hFILE *fpr;
        if ((fpr = hdopen(fd, mode)) == 0) return 0;
        fp = bgzf_read_init(fpr);
        if (fp == 0) { hclose_abruptly(fpr); return NULL; } // FIXME this closes fd
        fp->fp = fpr;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        hFILE *fpw;
        if ((fpw = hdopen(fd, mode)) == 0) return 0;
        fp = bgzf_write_init(mode);
        if (fp == NULL) return NULL;
        fp->fp = fpw;
    }
    else { errno = EINVAL; return 0; }

    fp->is_be = ed_is_big();
    return fp;
}

BGZF *bgzf_hopen(hFILE *hfp, const char *mode)
{
    BGZF *fp = NULL;
    assert(compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);
    if (strchr(mode, 'r')) {
        fp = bgzf_read_init(hfp);
        if (fp == NULL) return NULL;
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        fp = bgzf_write_init(mode);
        if (fp == NULL) return NULL;
    }
    else { errno = EINVAL; return 0; }

    fp->fp = hfp;
    fp->is_be = ed_is_big();
    return fp;
}

int bgzf_compress(void *_dst, size_t *dlen, const void *src, size_t slen, int level)
{
    uint32_t crc;
    z_stream zs;
    uint8_t *dst = (uint8_t*)_dst;

    // compress the body
    zs.zalloc = NULL; zs.zfree = NULL;
    zs.msg = NULL;
    zs.next_in  = (Bytef*)src;
    zs.avail_in = slen;
    zs.next_out = dst + BLOCK_HEADER_LENGTH;
    zs.avail_out = *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
    int ret = deflateInit2(&zs, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY); // -15 to disable zlib header/footer
    if (ret!=Z_OK) {
        hts_log_error("Call to deflateInit2 failed: %s", bgzf_zerr(ret, &zs));
        return -1;
    }
    if ((ret = deflate(&zs, Z_FINISH)) != Z_STREAM_END) {
        hts_log_error("Deflate operation failed: %s", bgzf_zerr(ret, ret == Z_DATA_ERROR ? &zs : NULL));
        return -1;
    }
    if ((ret = deflateEnd(&zs)) != Z_OK) {
        hts_log_error("Call to deflateEnd failed: %s", bgzf_zerr(ret, NULL));
        return -1;
    }
    *dlen = zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
    // write the header
    memcpy(dst, g_magic, BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
    packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes
    // write the footer
    crc = crc32(crc32(0L, NULL, 0L), (Bytef*)src, slen);
    packInt32((uint8_t*)&dst[*dlen - 8], crc);
    packInt32((uint8_t*)&dst[*dlen - 4], slen);
    return 0;
}

static int bgzf_gzip_compress(BGZF *fp, void *_dst, size_t *dlen, const void *src, size_t slen, int level)
{
    uint8_t *dst = (uint8_t*)_dst;
    z_stream *zs = fp->gz_stream;
    int flush = slen ? Z_PARTIAL_FLUSH : Z_FINISH;
    zs->next_in   = (Bytef*)src;
    zs->avail_in  = slen;
    zs->next_out  = dst;
    zs->avail_out = *dlen;
    int ret = deflate(zs, flush);
    if (ret == Z_STREAM_ERROR) {
        hts_log_error("Deflate operation failed: %s", bgzf_zerr(ret, NULL));
        return -1;
    }
    if (zs->avail_in != 0) {
        hts_log_error("Deflate block too large for output buffer");
        return -1;
    }
    *dlen = *dlen - zs->avail_out;
    return 0;
}

// Deflate the block in fp->uncompressed_block into fp->compressed_block. Also adds an extra field that stores the compressed block length.
static int deflate_block(BGZF *fp, int block_length)
{
    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret;
    if ( !fp->is_gzip )
        ret = bgzf_compress(fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level);
    else
        ret = bgzf_gzip_compress(fp, fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level);

    if ( ret != 0 )
    {
        hts_log_debug("Compression error %d", ret);
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    fp->block_offset = 0;
    return comp_size;
}

static int bgzf_uncompress(uint8_t *dst, size_t *dlen, const uint8_t *src, size_t slen) {
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.msg = NULL;
    zs.next_in = (Bytef*)src;
    zs.avail_in = slen;
    zs.next_out = (Bytef*)dst;
    zs.avail_out = *dlen;

    int ret = inflateInit2(&zs, -15);
    if (ret != Z_OK) {
        hts_log_error("Call to inflateInit2 failed: %s", bgzf_zerr(ret, &zs));
        return -1;
    }
    if ((ret = inflate(&zs, Z_FINISH)) != Z_STREAM_END) {
        hts_log_error("Inflate operation failed: %s", bgzf_zerr(ret, ret == Z_DATA_ERROR ? &zs : NULL));
        if ((ret = inflateEnd(&zs)) != Z_OK) {
            hts_log_warning("Call to inflateEnd failed: %s", bgzf_zerr(ret, NULL));
        }
        return -1;
    }
    if ((ret = inflateEnd(&zs)) != Z_OK) {
        hts_log_error("Call to inflateEnd failed: %s", bgzf_zerr(ret, NULL));
        return -1;
    }
    *dlen = *dlen - zs.avail_out;
    return 0;
}

// Inflate the block in fp->compressed_block into fp->uncompressed_block
static int inflate_block(BGZF* fp, int block_length)
{
    size_t dlen = BGZF_MAX_BLOCK_SIZE;
    int ret = bgzf_uncompress(fp->uncompressed_block, &dlen,
                              (Bytef*)fp->compressed_block + 18, block_length - 18);
    if (ret < 0) {
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }

    // Check CRC of uncompressed block matches the gzip header.
    // NB: we may wish to switch out the zlib crc32 for something more performant.
    // See PR#361 and issue#467
    uint32_t c1 = crc32(0L, (unsigned char *)fp->uncompressed_block, dlen);
    uint32_t c2 = le_to_u32((uint8_t *)fp->compressed_block + block_length-8);
    if (c1 != c2) {
        fp->errcode |= BGZF_ERR_CRC;
        return -1;
    }

    return dlen;
}

static int inflate_gzip_block(BGZF *fp, int cached)
{
    int ret = Z_OK;
    do
    {
        if ( !cached && fp->gz_stream->avail_out!=0 )
        {
            fp->gz_stream->avail_in = hread(fp->fp, fp->compressed_block, BGZF_BLOCK_SIZE);
            if ( fp->gz_stream->avail_in<=0 ) return fp->gz_stream->avail_in;
            if ( fp->gz_stream->avail_in==0 ) break;
            fp->gz_stream->next_in = fp->compressed_block;
        }
        else cached = 0;
        do
        {
            fp->gz_stream->next_out = (Bytef*)fp->uncompressed_block + fp->block_offset;
            fp->gz_stream->avail_out = BGZF_MAX_BLOCK_SIZE - fp->block_offset;
            fp->gz_stream->msg = NULL;
            ret = inflate(fp->gz_stream, Z_NO_FLUSH);
            if (ret < 0 && ret != Z_BUF_ERROR) {
                hts_log_error("Inflate operation failed: %s", bgzf_zerr(ret, ret == Z_DATA_ERROR ? fp->gz_stream : NULL));
                fp->errcode |= BGZF_ERR_ZLIB;
                return -1;
            }
            unsigned int have = BGZF_MAX_BLOCK_SIZE - fp->gz_stream->avail_out;
            if ( have ) return have;
        }
        while ( fp->gz_stream->avail_out == 0 );
    }
    while (ret != Z_STREAM_END);
    return BGZF_MAX_BLOCK_SIZE - fp->gz_stream->avail_out;
}

// Returns: 0 on success (BGZF header); -1 on non-BGZF GZIP header; -2 on error
static int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}

#ifdef BGZF_CACHE
static void free_cache(BGZF *fp)
{
    khint_t k;
    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    if (fp->is_write) return;
    for (k = kh_begin(h); k < kh_end(h); ++k)
        if (kh_exist(h, k)) free(kh_val(h, k).block);
    kh_destroy(cache, h);
}

static int load_block_from_cache(BGZF *fp, int64_t block_address)
{
    khint_t k;
    cache_t *p;

    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    k = kh_get(cache, h, block_address);
    if (k == kh_end(h)) return 0;
    p = &kh_val(h, k);
    if (fp->block_length != 0) fp->block_offset = 0;
    fp->block_address = block_address;
    fp->block_length = p->size;
    // FIXME: why BGZF_MAX_BLOCK_SIZE and not p->size?
    memcpy(fp->uncompressed_block, p->block, BGZF_MAX_BLOCK_SIZE);
    if ( hseek(fp->fp, p->end_offset, SEEK_SET) < 0 )
    {
        // todo: move the error up
        hts_log_error("Could not hseek to %"PRId64"", p->end_offset);
        exit(1);
    }
    return p->size;
}

static void cache_block(BGZF *fp, int size)
{
    int ret;
    khint_t k;
    cache_t *p;
    //fprintf(stderr, "Cache block at %llx\n", (int)fp->block_address);
    khash_t(cache) *h = (khash_t(cache)*)fp->cache;
    if (BGZF_MAX_BLOCK_SIZE >= fp->cache_size) return;
    if ((kh_size(h) + 1) * BGZF_MAX_BLOCK_SIZE > (uint32_t)fp->cache_size) {
        /* A better way would be to remove the oldest block in the
         * cache, but here we remove a random one for simplicity. This
         * should not have a big impact on performance. */
        for (k = kh_begin(h); k < kh_end(h); ++k)
            if (kh_exist(h, k)) break;
        if (k < kh_end(h)) {
            free(kh_val(h, k).block);
            kh_del(cache, h, k);
        }
    }
    k = kh_put(cache, h, fp->block_address, &ret);
    if (ret == 0) return; // if this happens, a bug!
    p = &kh_val(h, k);
    p->size = fp->block_length;
    p->end_offset = fp->block_address + size;
    p->block = (uint8_t*)malloc(BGZF_MAX_BLOCK_SIZE);
    memcpy(kh_val(h, k).block, fp->uncompressed_block, BGZF_MAX_BLOCK_SIZE);
}
#else
static void free_cache(BGZF *fp) {}
static int load_block_from_cache(BGZF *fp, int64_t block_address) {return 0;}
static void cache_block(BGZF *fp, int size) {}
#endif

/*
 * Absolute htell in this compressed file.
 *
 * Do not confuse with the external bgzf_tell macro which returns the virtual
 * offset.
 */
static off_t bgzf_htell(BGZF *fp) {
    if (fp->mt) {
        pthread_mutex_lock(&fp->mt->job_pool_m);
        off_t pos = fp->block_address + fp->block_clength;
        pthread_mutex_unlock(&fp->mt->job_pool_m);
        return pos;
    } else {
        return htell(fp->fp);
    }
}

int bgzf_read_block(BGZF *fp)
{
    hts_tpool_result *r;

    if (fp->mt) {
    again:
        if (fp->mt->hit_eof) {
            // Further reading at EOF will always return 0
            fp->block_length = 0;
            return 0;
        }
        r = hts_tpool_next_result_wait(fp->mt->out_queue);
        bgzf_job *j = r ? (bgzf_job *)hts_tpool_result_data(r) : NULL;

        if (!j || j->errcode == BGZF_ERR_MT) {
            if (!fp->mt->free_block) {
                fp->uncompressed_block = malloc(2 * BGZF_MAX_BLOCK_SIZE);
                if (fp->uncompressed_block == NULL) return -1;
                fp->compressed_block = (char *)fp->uncompressed_block + BGZF_MAX_BLOCK_SIZE;
            } // else it's already allocated with malloc, maybe even in-use.
            mt_destroy(fp->mt);
            fp->mt = NULL;
            hts_tpool_delete_result(r, 0);

            goto single_threaded;
        }

        if (j->errcode) {
            fp->errcode = j->errcode;
            return -1;
        }

        if (j->hit_eof) {
            if (!fp->last_block_eof && !fp->no_eof_block) {
                fp->no_eof_block = 1;
                hts_log_warning("EOF marker is absent. The input is probably truncated");
            }
            fp->mt->hit_eof = 1;
        }

        // Zero length blocks in the middle of a file are (wrongly)
        // considered as EOF by many callers.  We work around this by
        // trying again to see if we hit a genuine EOF.
        if (!j->hit_eof && j->uncomp_len == 0) {
            fp->last_block_eof = 1;
            hts_tpool_delete_result(r, 0);
            goto again;
        }

        // block_length=0 and block_offset set by bgzf_seek.
        if (fp->block_length != 0) fp->block_offset = 0;
        if (!j->hit_eof) fp->block_address = j->block_address;
        fp->block_clength = j->comp_len;
        fp->block_length = j->uncomp_len;
        // bgzf_read() can change fp->block_length
        fp->last_block_eof = (fp->block_length == 0);

        if ( j->uncomp_len && j->fp->idx_build_otf )
        {
            bgzf_index_add_block(j->fp);
            j->fp->idx->ublock_addr += j->uncomp_len;
        }

        // Steal the data block as it's quicker than a memcpy.
        // We just need to make sure we delay the pool free.
        if (fp->mt->curr_job) {
            pthread_mutex_lock(&fp->mt->job_pool_m);
            pool_free(fp->mt->job_pool, fp->mt->curr_job);
            pthread_mutex_unlock(&fp->mt->job_pool_m);
        }
        fp->uncompressed_block = j->uncomp_data;
        fp->mt->curr_job = j;
        if (fp->mt->free_block) {
            free(fp->mt->free_block); // clear up last non-mt block
            fp->mt->free_block = NULL;
        }

        hts_tpool_delete_result(r, 0);
        return 0;
    }

    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size, block_length, remaining;

 single_threaded:
    size = 0;

    // Reading an uncompressed file
    if ( !fp->is_compressed )
    {
        count = hread(fp->fp, fp->uncompressed_block, BGZF_MAX_BLOCK_SIZE);
        if (count < 0)  // Error
        {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        else if (count == 0)  // EOF
        {
            fp->block_length = 0;
            return 0;
        }
        if (fp->block_length != 0) fp->block_offset = 0;
        fp->block_address += count;
        fp->block_length = count;
        return 0;
    }

    // Reading compressed file
    int64_t block_address;
    block_address = bgzf_htell(fp);
    if ( fp->is_gzip && fp->gz_stream ) // is this is an initialized gzip stream?
    {
        count = inflate_gzip_block(fp, 0);
        if ( count<0 )
        {
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->block_length = count;
        fp->block_address = block_address;
        return 0;
    }
    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;

    // loop to skip empty bgzf blocks
    while (1)
    {
        count = hread(fp->fp, header, sizeof(header));
        if (count == 0) { // no data read
            if (!fp->last_block_eof && !fp->no_eof_block && !fp->is_gzip) {
                fp->no_eof_block = 1;
                hts_log_warning("EOF marker is absent. The input is probably truncated");
            }
            fp->block_length = 0;
            return 0;
        }
        int ret;
        if ( count != sizeof(header) || (ret=check_header(header))==-2 )
        {
            fp->errcode |= BGZF_ERR_HEADER;
            return -1;
        }
        if ( ret==-1 )
        {
            // GZIP, not BGZF
            uint8_t *cblock = (uint8_t*)fp->compressed_block;
            memcpy(cblock, header, sizeof(header));
            count = hread(fp->fp, cblock+sizeof(header), BGZF_BLOCK_SIZE - sizeof(header)) + sizeof(header);
            int nskip = 10;

            // Check optional fields to skip: FLG.FNAME,FLG.FCOMMENT,FLG.FHCRC,FLG.FEXTRA
            // Note: Some of these fields are untested, I did not have appropriate data available
            if ( header[3] & 0x4 ) // FLG.FEXTRA
            {
                nskip += unpackInt16(&cblock[nskip]) + 2;
            }
            if ( header[3] & 0x8 ) // FLG.FNAME
            {
                while ( nskip<count && cblock[nskip] ) nskip++;
                nskip++;
            }
            if ( header[3] & 0x10 ) // FLG.FCOMMENT
            {
                while ( nskip<count && cblock[nskip] ) nskip++;
                nskip++;
            }
            if ( header[3] & 0x2 ) nskip += 2;  //  FLG.FHCRC

            /* FIXME: Should handle this better.  There's no reason why
               someone shouldn't include a massively long comment in their
               gzip stream. */
            if ( nskip >= count )
            {
                fp->errcode |= BGZF_ERR_HEADER;
                return -1;
            }

            fp->is_gzip = 1;
            fp->gz_stream = (z_stream*) calloc(1,sizeof(z_stream));
            int ret = inflateInit2(fp->gz_stream, -15);
            if (ret != Z_OK)
            {
                hts_log_error("Call to inflateInit2 failed: %s", bgzf_zerr(ret, fp->gz_stream));
                fp->errcode |= BGZF_ERR_ZLIB;
                return -1;
            }
            fp->gz_stream->avail_in = count - nskip;
            fp->gz_stream->next_in  = cblock + nskip;
            count = inflate_gzip_block(fp, 1);
            if ( count<0 )
            {
                fp->errcode |= BGZF_ERR_ZLIB;
                return -1;
            }
            fp->block_length = count;
            fp->block_address = block_address;
            if ( fp->idx_build_otf ) return -1; // cannot build index for gzip
            return 0;
        }
        size = count;
        block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
        if (block_length < BLOCK_HEADER_LENGTH)
        {
            fp->errcode |= BGZF_ERR_HEADER;
            return -1;
        }
        compressed_block = (uint8_t*)fp->compressed_block;
        memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
        remaining = block_length - BLOCK_HEADER_LENGTH;
        count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
        if (count != remaining) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        size += count;
        if ((count = inflate_block(fp, block_length)) < 0) {
            hts_log_debug("Inflate block operation failed: %s", bgzf_zerr(count, NULL));
            fp->errcode |= BGZF_ERR_ZLIB;
            return -1;
        }
        fp->last_block_eof = (count == 0);
        if ( count ) break;     // otherwise an empty bgzf block
    }
    if (fp->block_length != 0) fp->block_offset = 0; // Do not reset offset if this read follows a seek.
    fp->block_address = block_address;
    fp->block_length = count;
    if ( fp->idx_build_otf )
    {
        bgzf_index_add_block(fp);
        fp->idx->ublock_addr += count;
    }
    cache_block(fp, size);
    return 0;
}

ssize_t bgzf_read(BGZF *fp, void *data, size_t length)
{
    ssize_t bytes_read = 0;
    uint8_t *output = (uint8_t*)data;
    if (length <= 0) return 0;
    assert(fp->is_write == 0);
    while (bytes_read < length) {
        int copy_length, available = fp->block_length - fp->block_offset;
        uint8_t *buffer;
        if (available <= 0) {
            int ret = bgzf_read_block(fp);
            if (ret != 0) {
                hts_log_error("Read block operation failed with error %d after %zd of %zu bytes", ret, bytes_read, length);
                fp->errcode |= BGZF_ERR_ZLIB;
                return -1;
            }
            available = fp->block_length - fp->block_offset;
            if (available <= 0) break;
        }
        copy_length = length - bytes_read < available? length - bytes_read : available;
        buffer = (uint8_t*)fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;

        // For raw gzip streams this avoids short reads.
        if (fp->block_offset == fp->block_length) {
            fp->block_address = bgzf_htell(fp);
            fp->block_offset = fp->block_length = 0;
        }
    }

    fp->uncompressed_address += bytes_read;

    return bytes_read;
}

ssize_t bgzf_raw_read(BGZF *fp, void *data, size_t length)
{
    ssize_t ret = hread(fp->fp, data, length);
    if (ret < 0) fp->errcode |= BGZF_ERR_IO;
    return ret;
}

#ifdef BGZF_MT

void *bgzf_encode_func(void *arg) {
    bgzf_job *j = (bgzf_job *)arg;

    j->comp_len = BGZF_MAX_BLOCK_SIZE;
    int ret = bgzf_compress(j->comp_data, &j->comp_len,
                            j->uncomp_data, j->uncomp_len,
                            j->fp->compress_level);
    if (ret != 0)
        j->errcode |= BGZF_ERR_ZLIB;

    return arg;
}

// Our input block has already been decoded by bgzf_mt_read_block().
// We need to split that into a fetch block (compressed) and make this
// do the actual decompression step.
void *bgzf_decode_func(void *arg) {
    bgzf_job *j = (bgzf_job *)arg;

    j->uncomp_len = BGZF_MAX_BLOCK_SIZE;
    int ret = bgzf_uncompress(j->uncomp_data, &j->uncomp_len,
                              j->comp_data+18, j->comp_len-18);
    if (ret != 0)
        j->errcode |= BGZF_ERR_ZLIB;

    return arg;
}

/*
 * Nul function so we can dispatch a job with the correct serial
 * to mark failure or to indicate an empty read (EOF).
 */
void *bgzf_nul_func(void *arg) { return arg; }

/*
 * Takes compressed blocks off the results queue and calls hwrite to
 * punt them to the output stream.
 *
 * Returns NULL when no more are left, or -1 on error
 */
static void *bgzf_mt_writer(void *vp) {
    BGZF *fp = (BGZF *)vp;
    mtaux_t *mt = fp->mt;
    hts_tpool_result *r;

    // Iterates until result queue is shutdown, where it returns NULL.
    while ((r = hts_tpool_next_result_wait(mt->out_queue))) {
        bgzf_job *j = (bgzf_job *)hts_tpool_result_data(r);
        assert(j);

        if (hwrite(fp->fp, j->comp_data, j->comp_len) != j->comp_len) {
            fp->errcode |= BGZF_ERR_IO;
            goto err;
        }

        /*
         * Periodically call hflush (which calls fsync when on a file).
         * This avoids the fsync being done at the bgzf_close stage,
         * which can sometimes cause signficant delays.  As this is in
         * a separate thread, spreading the sync delays throughout the
         * program execution seems better.
         * Frequency of 1/512 has been chosen by experimentation
         * across local XFS, NFS and Lustre tests.
         */
        if (++mt->flush_pending % 512 == 0)
            if (hflush(fp->fp) != 0)
                goto err;


        hts_tpool_delete_result(r, 0);

        // Also updated by main thread
        pthread_mutex_lock(&mt->job_pool_m);
        pool_free(mt->job_pool, j);
        mt->jobs_pending--;
        pthread_mutex_unlock(&mt->job_pool_m);
    }

    if (hflush(fp->fp) != 0)
        goto err;

    hts_tpool_process_destroy(mt->out_queue);

    return NULL;

 err:
    hts_tpool_process_destroy(mt->out_queue);
    return (void *)-1;
}


/*
 * Reads a compressed block of data using hread and dispatches it to
 * the thread pool for decompression.  This is the analogue of the old
 * non-threaded bgzf_read_block() function, but without modifying fp
 * in any way (except for the read offset).  All output goes via the
 * supplied bgzf_job struct.
 *
 * Returns NULL when no more are left, or -1 on error
 */
int bgzf_mt_read_block(BGZF *fp, bgzf_job *j)
{
    uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
    int count, size = 0, block_length, remaining;

    // NOTE: Guaranteed to be compressed as we block multi-threading in
    // uncompressed mode.  However it may be gzip compression instead
    // of bgzf.

    // Reading compressed file
    int64_t block_address;
    block_address = htell(fp->fp);

    if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
    count = hpeek(fp->fp, header, sizeof(header));
    if (count == 0) // no data read
        return -1;
    int ret;
    if ( count != sizeof(header) || (ret=check_header(header))==-2 )
    {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    if (ret == -1) {
        j->errcode |= BGZF_ERR_MT;
        return -1;
    }

    count = hread(fp->fp, header, sizeof(header));
    if (count != sizeof(header)) // no data read
        return -1;

    size = count;
    block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
    if (block_length < BLOCK_HEADER_LENGTH) {
        j->errcode |= BGZF_ERR_HEADER;
        return -1;
    }
    compressed_block = (uint8_t*)j->comp_data;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    remaining = block_length - BLOCK_HEADER_LENGTH;
    count = hread(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
    if (count != remaining) {
        j->errcode |= BGZF_ERR_IO;
        return -1;
    }
    size += count;
    j->comp_len = block_length;
    j->uncomp_len = BGZF_MAX_BLOCK_SIZE;
    j->block_address = block_address;
    j->fp = fp;
    j->errcode = 0;

    return 0;
}


static int bgzf_check_EOF_common(BGZF *fp)
{
    uint8_t buf[28];
    off_t offset = htell(fp->fp);
    if (hseek(fp->fp, -28, SEEK_END) < 0) {
        if (errno == ESPIPE) { hclearerr(fp->fp); return 2; }
#ifdef _WIN32
        if (errno == EINVAL) { hclearerr(fp->fp); return 2; }
#endif
        else return -1;
    }
    if ( hread(fp->fp, buf, 28) != 28 ) return -1;
    if ( hseek(fp->fp, offset, SEEK_SET) < 0 ) return -1;
    return (memcmp("\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0", buf, 28) == 0)? 1 : 0;
}

/*
 * Checks EOF from the reader thread.
 */
static void bgzf_mt_eof(BGZF *fp) {
    mtaux_t *mt = fp->mt;

    pthread_mutex_lock(&mt->job_pool_m);
    mt->eof = bgzf_check_EOF_common(fp);
    pthread_mutex_unlock(&mt->job_pool_m);
    pthread_cond_signal(&mt->command_c);
}


/*
 * Performs the seek (called by reader thread).
 *
 * This simply drains the entire queue, throwing away blocks, seeks,
 * and starts it up again.  Brute force, but maybe sufficient.
 */
static void bgzf_mt_seek(BGZF *fp) {
    mtaux_t *mt = fp->mt;

    hts_tpool_process_reset(mt->out_queue, 0);
    pthread_mutex_lock(&mt->job_pool_m);
    mt->command = NONE;
    mt->errcode = 0;

    if (hseek(fp->fp, mt->block_address, SEEK_SET) < 0)
        mt->errcode = BGZF_ERR_IO;

    pthread_mutex_unlock(&mt->job_pool_m);
    pthread_cond_signal(&mt->command_c);
}

static void *bgzf_mt_reader(void *vp) {
    BGZF *fp = (BGZF *)vp;
    mtaux_t *mt = fp->mt;

restart:
    pthread_mutex_lock(&mt->job_pool_m);
    bgzf_job *j = pool_alloc(mt->job_pool);
    pthread_mutex_unlock(&mt->job_pool_m);
    j->errcode = 0;
    j->comp_len = 0;
    j->uncomp_len = 0;
    j->hit_eof = 0;

    while (bgzf_mt_read_block(fp, j) == 0) {
        // Dispatch
        hts_tpool_dispatch(mt->pool, mt->out_queue, bgzf_decode_func, j);

        // Check for command
        pthread_mutex_lock(&mt->command_m);
        switch (mt->command) {
        case SEEK:
            bgzf_mt_seek(fp);
            pthread_mutex_unlock(&mt->command_m);
            goto restart;

        case HAS_EOF:
            bgzf_mt_eof(fp);
            break;

        case CLOSE:
            pthread_cond_signal(&mt->command_c);
            pthread_mutex_unlock(&mt->command_m);
            hts_tpool_process_destroy(mt->out_queue);
            pthread_exit(NULL);

        default:
            break;
        }
        pthread_mutex_unlock(&mt->command_m);

        // Allocate buffer for next block
        pthread_mutex_lock(&mt->job_pool_m);
        j = pool_alloc(mt->job_pool);
        pthread_mutex_unlock(&mt->job_pool_m);
        j->errcode = 0;
        j->comp_len = 0;
        j->uncomp_len = 0;
        j->hit_eof = 0;
    }

    if (j->errcode == BGZF_ERR_MT) {
        // Attempt to multi-thread decode a raw gzip stream cannot be done.
        // We tear down the multi-threaded decoder and revert to the old code.
        hts_tpool_dispatch(mt->pool, mt->out_queue, bgzf_nul_func, j);
        hts_tpool_process_ref_decr(mt->out_queue);
        pthread_exit(&j->errcode);
    }

    // Dispatch an empty block so EOF is spotted.
    // We also use this mechanism for returning errors, in which case
    // j->errcode is set already.

    j->hit_eof = 1;
    hts_tpool_dispatch(mt->pool, mt->out_queue, bgzf_nul_func, j);
    if (j->errcode != 0) {
        hts_tpool_process_destroy(mt->out_queue);
        pthread_exit(&j->errcode);
    }

    // We hit EOF so can stop reading, but we may get a subsequent
    // seek request.  In this case we need to restart the reader.
    //
    // To handle this we wait on a condition variable and then
    // monitor the command. (This could be either seek or close.)
    for (;;) {
        pthread_mutex_lock(&mt->command_m);
        if (mt->command == NONE)
            pthread_cond_wait(&mt->command_c, &mt->command_m);
        switch(mt->command) {
        default:
            pthread_mutex_unlock(&mt->command_m);
            break;

        case SEEK:
            bgzf_mt_seek(fp);
            pthread_mutex_unlock(&mt->command_m);
            goto restart;

        case HAS_EOF:
            bgzf_mt_eof(fp);
            pthread_mutex_unlock(&mt->command_m);
            continue;

        case CLOSE:
            pthread_cond_signal(&mt->command_c);
            pthread_mutex_unlock(&mt->command_m);
            hts_tpool_process_destroy(mt->out_queue);
            pthread_exit(NULL);
        }
    }
    return NULL;
}

int bgzf_thread_pool(BGZF *fp, hts_tpool *pool, int qsize) {
    // No gain from multi-threading when not compressed
    if (!fp->is_compressed)
        return 0;

    mtaux_t *mt;
    mt = (mtaux_t*)calloc(1, sizeof(mtaux_t));
    if (!mt) return -1;
    fp->mt = mt;

    mt->pool = pool;
    mt->n_threads = hts_tpool_size(pool);
    if (!qsize)
        qsize = mt->n_threads*2;
    if (!(mt->out_queue = hts_tpool_process_init(mt->pool, qsize, 0))) {
        free(mt);
        return -1;
    }
    hts_tpool_process_ref_incr(mt->out_queue);

    mt->job_pool = pool_create(sizeof(bgzf_job));

    pthread_mutex_init(&mt->job_pool_m, NULL);
    pthread_mutex_init(&mt->command_m, NULL);
    pthread_cond_init(&mt->command_c, NULL);
    mt->flush_pending = 0;
    mt->jobs_pending = 0;
    mt->free_block = fp->uncompressed_block; // currently in-use block
    pthread_create(&mt->io_task, NULL,
                   fp->is_write ? bgzf_mt_writer : bgzf_mt_reader, fp);

    return 0;
}

int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)
{
    // No gain from multi-threading when not compressed
    if (!fp->is_compressed || fp->is_gzip)
        return 0;

    if (n_threads < 1) return -1;
    hts_tpool *p = hts_tpool_init(n_threads);
    if (!p)
        return -1;

    if (bgzf_thread_pool(fp, p, 0) != 0) {
        hts_tpool_destroy(p);
        return -1;
    }

    fp->mt->own_pool = 1;

    return 0;
}

static void mt_destroy(mtaux_t *mt)
{
    pthread_mutex_lock(&mt->command_m);
    mt->command = CLOSE;
    pthread_cond_signal(&mt->command_c);
    hts_tpool_wake_dispatch(mt->out_queue); // unstick the reader
    pthread_mutex_unlock(&mt->command_m);

    // Destroying the queue first forces the writer to exit.
    hts_tpool_process_destroy(mt->out_queue);
    pthread_join(mt->io_task, NULL);

    pthread_mutex_destroy(&mt->job_pool_m);
    pthread_mutex_destroy(&mt->command_m);
    pthread_cond_destroy(&mt->command_c);
    if (mt->curr_job)
        pool_free(mt->job_pool, mt->curr_job);

    if (mt->own_pool)
        hts_tpool_destroy(mt->pool);

    pool_destroy(mt->job_pool);

    free(mt);
    fflush(stderr);
}

static int mt_queue(BGZF *fp)
{
    mtaux_t *mt = fp->mt;

    // Also updated by writer thread
    pthread_mutex_lock(&mt->job_pool_m);
    bgzf_job *j = pool_alloc(mt->job_pool);
    mt->jobs_pending++;
    pthread_mutex_unlock(&mt->job_pool_m);

    j->fp = fp;
    j->errcode = 0;
    j->uncomp_len  = fp->block_offset;
    memcpy(j->uncomp_data, fp->uncompressed_block, j->uncomp_len);

    // Need non-block vers & job_pending?
    hts_tpool_dispatch(mt->pool, mt->out_queue, bgzf_encode_func, j);

    fp->block_offset = 0;
    return 0;
}

static int mt_flush_queue(BGZF *fp)
{
    mtaux_t *mt = fp->mt;

    // Drain the encoder jobs.
    // We cannot use hts_tpool_flush here as it can cause deadlock if
    // the queue is full up of decoder tasks.  The best solution would
    // be to have one input queue per type of job, but we don't right now.
    //hts_tpool_flush(mt->pool);
    pthread_mutex_lock(&mt->job_pool_m);
    while (mt->jobs_pending != 0) {
        pthread_mutex_unlock(&mt->job_pool_m);
        usleep(10000); // FIXME: replace by condition variable
        pthread_mutex_lock(&mt->job_pool_m);
    }
    pthread_mutex_unlock(&mt->job_pool_m);

    // Wait on bgzf_mt_writer to drain the queue
    if (hts_tpool_process_flush(mt->out_queue) != 0)
        return -1;

    return (fp->errcode == 0)? 0 : -1;
}

static int lazy_flush(BGZF *fp)
{
    if (fp->mt)
        return fp->block_offset ? mt_queue(fp) : 0;
    else
        return bgzf_flush(fp);
}

#else  // ~ #ifdef BGZF_MT

int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)
{
    return 0;
}

static inline int lazy_flush(BGZF *fp)
{
    return bgzf_flush(fp);
}

#endif // ~ #ifdef BGZF_MT

int bgzf_flush(BGZF *fp)
{
    if (!fp->is_write) return 0;
#ifdef BGZF_MT
    if (fp->mt) {
        int ret = 0;
        if (fp->block_offset) ret = mt_queue(fp);
        return ret ? ret : mt_flush_queue(fp);
    }
#endif
    while (fp->block_offset > 0) {
        int block_length;
        if ( fp->idx_build_otf )
        {
            bgzf_index_add_block(fp);
            fp->idx->ublock_addr += fp->block_offset;
        }
        block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, fp->compressed_block, block_length) != block_length) {
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }
        fp->block_address += block_length;
    }
    return 0;
}

int bgzf_flush_try(BGZF *fp, ssize_t size)
{
    if (fp->block_offset + size > BGZF_BLOCK_SIZE) return lazy_flush(fp);
    return 0;
}

ssize_t bgzf_write(BGZF *fp, const void *data, size_t length)
{
    if ( !fp->is_compressed )
        return hwrite(fp->fp, data, length);

    const uint8_t *input = (const uint8_t*)data;
    ssize_t remaining = length;
    assert(fp->is_write);
    while (remaining > 0) {
        uint8_t* buffer = (uint8_t*)fp->uncompressed_block;
        int copy_length = BGZF_BLOCK_SIZE - fp->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (fp->block_offset == BGZF_BLOCK_SIZE) {
            if (lazy_flush(fp) != 0) return -1;
        }
    }
    return length - remaining;
}

ssize_t bgzf_block_write(BGZF *fp, const void *data, size_t length)
{
    if ( !fp->is_compressed )
        return hwrite(fp->fp, data, length);
    const uint8_t *input = (const uint8_t*)data;
    ssize_t remaining = length;
    assert(fp->is_write);
    uint64_t current_block; //keep track of current block
    uint64_t ublock_size; // amount of uncompressed data to be fed into next block
    while (remaining > 0) {
        current_block = fp->idx->moffs - fp->idx->noffs;
        ublock_size = fp->idx->offs[current_block+1].uaddr-fp->idx->offs[current_block].uaddr;
        uint8_t* buffer = (uint8_t*)fp->uncompressed_block;
        int copy_length = ublock_size - fp->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (fp->block_offset == ublock_size) {
            if (lazy_flush(fp) != 0) return -1;
            fp->idx->noffs--;  // decrement noffs to track the blocks
        }
    }
    return length - remaining;
}


ssize_t bgzf_raw_write(BGZF *fp, const void *data, size_t length)
{
    ssize_t ret = hwrite(fp->fp, data, length);
    if (ret < 0) fp->errcode |= BGZF_ERR_IO;
    return ret;
}

int bgzf_close(BGZF* fp)
{
    int ret, block_length;
    if (fp == 0) return -1;
    if (fp->is_write && fp->is_compressed) {
        if (bgzf_flush(fp) != 0) return -1;
        fp->compress_level = -1;
        block_length = deflate_block(fp, 0); // write an empty block
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, fp->compressed_block, block_length) < 0
            || hflush(fp->fp) != 0) {
            hts_log_error("File write failed");
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
    }
#ifdef BGZF_MT
    if (fp->mt) {
        if (!fp->mt->free_block)
            fp->uncompressed_block = NULL;
        mt_destroy(fp->mt);
    }
#endif
    if ( fp->is_gzip )
    {
        if (fp->gz_stream == NULL) ret = Z_OK;
        else if (!fp->is_write) ret = inflateEnd(fp->gz_stream);
        else ret = deflateEnd(fp->gz_stream);
        if (ret != Z_OK) {
            hts_log_error("Call to inflateEnd/deflateEnd failed: %s", bgzf_zerr(ret, NULL));
        }
        free(fp->gz_stream);
    }
    ret = hclose(fp->fp);
    if (ret != 0) return -1;
    bgzf_index_destroy(fp);
    free(fp->uncompressed_block);
    free_cache(fp);
    free(fp);
    return 0;
}

void bgzf_set_cache_size(BGZF *fp, int cache_size)
{
    if (fp) fp->cache_size = cache_size;
}

int bgzf_check_EOF(BGZF *fp) {
    int has_eof;

    if (fp->mt) {
        pthread_mutex_lock(&fp->mt->command_m);
        fp->mt->command = HAS_EOF;
        pthread_cond_signal(&fp->mt->command_c);
        hts_tpool_wake_dispatch(fp->mt->out_queue);
        pthread_cond_wait(&fp->mt->command_c, &fp->mt->command_m);
        has_eof = fp->mt->eof;
        pthread_mutex_unlock(&fp->mt->command_m);
    } else {
        has_eof = bgzf_check_EOF_common(fp);
    }

    fp->no_eof_block = (has_eof == 0);

    return has_eof;
}

int64_t bgzf_seek(BGZF* fp, int64_t pos, int where)
{
    int block_offset;
    int64_t block_address;

    if (fp->is_write || where != SEEK_SET || fp->is_gzip) {
        fp->errcode |= BGZF_ERR_MISUSE;
        return -1;
    }
    block_offset = pos & 0xFFFF;
    block_address = pos >> 16;

    if (fp->mt) {
        // The reader runs asynchronous and does loops of:
        //    Read block
        //    Check & process command
        //    Dispatch decode job
        //
        // Once at EOF it then switches to loops of
        //    Wait for command
        //    Process command (possibly switching back to above loop).
        //
        // To seek we therefore send the reader thread a SEEK command,
        // waking it up if blocked in dispatch and signalling if
        // waiting for a command.  We then wait for the response so we
        // know the seek succeeded.
        pthread_mutex_lock(&fp->mt->command_m);
        fp->mt->hit_eof = 0;
        fp->mt->command = SEEK;
        fp->mt->block_address = block_address;
        pthread_cond_signal(&fp->mt->command_c);
        hts_tpool_wake_dispatch(fp->mt->out_queue);
        pthread_cond_wait(&fp->mt->command_c, &fp->mt->command_m);

        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = block_address;
        fp->block_offset = block_offset;

        pthread_mutex_unlock(&fp->mt->command_m);
    } else {
        if (hseek(fp->fp, block_address, SEEK_SET) < 0) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = block_address;
        fp->block_offset = block_offset;
    }

    return 0;
}

int bgzf_is_bgzf(const char *fn)
{
    uint8_t buf[16];
    int n;
    hFILE *fp;
    if ((fp = hopen(fn, "r")) == 0) return 0;
    n = hread(fp, buf, 16);
    if (hclose(fp) < 0) return 0;
    if (n != 16) return 0;
    return check_header(buf) == 0? 1 : 0;
}

int bgzf_compression(BGZF *fp)
{
    return (!fp->is_compressed)? no_compression : (fp->is_gzip)? gzip : bgzf;
}

int bgzf_getc(BGZF *fp)
{
    if (fp->block_offset+1 < fp->block_length) {
        fp->uncompressed_address++;
        return ((unsigned char*)fp->uncompressed_block)[fp->block_offset++];
    }

    int c;
    if (fp->block_offset >= fp->block_length) {
        if (bgzf_read_block(fp) != 0) return -2; /* error */
        if (fp->block_length == 0) return -1; /* end-of-file */
    }
    c = ((unsigned char*)fp->uncompressed_block)[fp->block_offset++];
    if (fp->block_offset == fp->block_length) {
        fp->block_address = bgzf_htell(fp);
        fp->block_offset = 0;
        fp->block_length = 0;
    }
    fp->uncompressed_address++;
    return c;
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

int bgzf_getline(BGZF *fp, int delim, kstring_t *str)
{
    int l, state = 0;
    str->l = 0;
    do {
        if (fp->block_offset >= fp->block_length) {
            if (bgzf_read_block(fp) != 0) { state = -2; break; }
            if (fp->block_length == 0) { state = -1; break; }
        }
        unsigned char *buf = fp->uncompressed_block;
        for (l = fp->block_offset; l < fp->block_length && buf[l] != delim; ++l);
        if (l < fp->block_length) state = 1;
        l -= fp->block_offset;
        if (str->l + l + 1 >= str->m) {
            str->m = str->l + l + 2;
            kroundup32(str->m);
            str->s = (char*)realloc(str->s, str->m);
        }
        memcpy(str->s + str->l, buf + fp->block_offset, l);
        str->l += l;
        fp->block_offset += l + 1;
        if (fp->block_offset >= fp->block_length) {
            fp->block_address = bgzf_htell(fp);
            fp->block_offset = 0;
            fp->block_length = 0;
        }
    } while (state == 0);
    if (str->l == 0 && state < 0) return state;
    fp->uncompressed_address += str->l + 1;
    if ( delim=='\n' && str->l>0 && str->s[str->l-1]=='\r' ) str->l--;
    str->s[str->l] = 0;
    return str->l;
}

void bgzf_index_destroy(BGZF *fp)
{
    if ( !fp->idx ) return;
    free(fp->idx->offs);
    free(fp->idx);
    fp->idx = NULL;
    fp->idx_build_otf = 0;
}

int bgzf_index_build_init(BGZF *fp)
{
    bgzf_index_destroy(fp);
    fp->idx = (bgzidx_t*) calloc(1,sizeof(bgzidx_t));
    if ( !fp->idx ) return -1;
    fp->idx_build_otf = 1;  // build index on the fly
    return 0;
}

int bgzf_index_add_block(BGZF *fp)
{
    fp->idx->noffs++;
    if ( fp->idx->noffs > fp->idx->moffs )
    {
        fp->idx->moffs = fp->idx->noffs;
        kroundup32(fp->idx->moffs);
        fp->idx->offs = (bgzidx1_t*) realloc(fp->idx->offs, fp->idx->moffs*sizeof(bgzidx1_t));
        if ( !fp->idx->offs ) return -1;
    }
    fp->idx->offs[ fp->idx->noffs-1 ].uaddr = fp->idx->ublock_addr;
    fp->idx->offs[ fp->idx->noffs-1 ].caddr = fp->block_address;
    return 0;
}

static inline int hwrite_uint64(uint64_t x, hFILE *f)
{
    if (ed_is_big()) x = ed_swap_8(x);
    if (hwrite(f, &x, sizeof(x)) != sizeof(x)) return -1;
    return 0;
}

static char * get_name_suffix(const char *bname, const char *suffix)
{
    size_t len = strlen(bname) + strlen(suffix) + 1;
    char *buff = malloc(len);
    if (!buff) return NULL;
    snprintf(buff, len, "%s%s", bname, suffix);
    return buff;
}

int bgzf_index_dump_hfile(BGZF *fp, struct hFILE *idx, const char *name)
{
    // Note that the index contains one extra record when indexing files opened
    // for reading. The terminating record is not present when opened for writing.
    // This is not a bug.

    int i;

    if (!fp->idx) {
        hts_log_error("Called for BGZF handle with no index");
        errno = EINVAL;
        return -1;
    }

    if (bgzf_flush(fp) != 0) return -1;

    if (hwrite_uint64(fp->idx->noffs - 1, idx) < 0) goto fail;
    for (i=1; i<fp->idx->noffs; i++)
    {
        if (hwrite_uint64(fp->idx->offs[i].caddr, idx) < 0) goto fail;
        if (hwrite_uint64(fp->idx->offs[i].uaddr, idx) < 0) goto fail;
    }
    return 0;

 fail:
    hts_log_error("Error writing to %s : %s", name ? name : "index", strerror(errno));
    return -1;
}

int bgzf_index_dump(BGZF *fp, const char *bname, const char *suffix)
{
    const char *name = bname, *msg = NULL;
    char *tmp = NULL;
    hFILE *idx = NULL;

    if (!fp->idx) {
        hts_log_error("Called for BGZF handle with no index");
        errno = EINVAL;
        return -1;
    }

    if ( suffix )
    {
        tmp = get_name_suffix(bname, suffix);
        if ( !tmp ) return -1;
        name = tmp;
    }

    idx = hopen(name, "wb");
    if ( !idx ) {
        msg = "Error opening";
        goto fail;
    }

    if (bgzf_index_dump_hfile(fp, idx, name) != 0) goto fail;

    if (hclose(idx) < 0)
    {
        idx = NULL;
        msg = "Error on closing";
        goto fail;
    }

    free(tmp);
    return 0;

 fail:
    if (msg != NULL) {
        hts_log_error("%s %s : %s", msg, name, strerror(errno));
    }
    if (idx) hclose_abruptly(idx);
    free(tmp);
    return -1;
}

static inline int hread_uint64(uint64_t *xptr, hFILE *f)
{
    if (hread(f, xptr, sizeof(*xptr)) != sizeof(*xptr)) return -1;
    if (ed_is_big()) ed_swap_8p(xptr);
    return 0;
}

int bgzf_index_load_hfile(BGZF *fp, struct hFILE *idx, const char *name)
{
    fp->idx = (bgzidx_t*) calloc(1,sizeof(bgzidx_t));
    if (fp->idx == NULL) goto fail;
    uint64_t x;
    if (hread_uint64(&x, idx) < 0) goto fail;

    fp->idx->noffs = fp->idx->moffs = x + 1;
    fp->idx->offs  = (bgzidx1_t*) malloc(fp->idx->moffs*sizeof(bgzidx1_t));
    if (fp->idx->offs == NULL) goto fail;
    fp->idx->offs[0].caddr = fp->idx->offs[0].uaddr = 0;

    int i;
    for (i=1; i<fp->idx->noffs; i++)
    {
        if (hread_uint64(&fp->idx->offs[i].caddr, idx) < 0) goto fail;
        if (hread_uint64(&fp->idx->offs[i].uaddr, idx) < 0) goto fail;
    }

    return 0;

 fail:
    hts_log_error("Error reading %s : %s", name ? name : "index", strerror(errno));
    if (fp->idx) {
        free(fp->idx->offs);
        free(fp->idx);
        fp->idx = NULL;
    }
    return -1;
}

int bgzf_index_load(BGZF *fp, const char *bname, const char *suffix)
{
    const char *name = bname, *msg = NULL;
    char *tmp = NULL;
    hFILE *idx = NULL;
    if ( suffix )
    {
        tmp = get_name_suffix(bname, suffix);
        if ( !tmp ) return -1;
        name = tmp;
    }

    idx = hopen(name, "rb");
    if ( !idx ) {
        msg = "Error opening";
        goto fail;
    }

    if (bgzf_index_load_hfile(fp, idx, name) != 0) goto fail;

    if (hclose(idx) != 0) {
        idx = NULL;
        msg = "Error closing";
        goto fail;
    }

    free(tmp);
    return 0;

 fail:
    if (msg != NULL) {
        hts_log_error("%s %s : %s", msg, name, strerror(errno));
    }
    if (idx) hclose_abruptly(idx);
    free(tmp);
    return -1;
}

int bgzf_useek(BGZF *fp, long uoffset, int where)
{
    if ( !fp->is_compressed )
    {
        if (hseek(fp->fp, uoffset, SEEK_SET) < 0)
        {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        fp->block_length = 0;  // indicates current block has not been loaded
        fp->block_address = uoffset;
        fp->block_offset = 0;
        if (bgzf_read_block(fp) < 0) {
            fp->errcode |= BGZF_ERR_IO;
            return -1;
        }
        fp->uncompressed_address = uoffset;
        return 0;
    }

    if ( !fp->idx )
    {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }

    // binary search
    int ilo = 0, ihi = fp->idx->noffs - 1;
    while ( ilo<=ihi )
    {
        int i = (ilo+ihi)*0.5;
        if ( uoffset < fp->idx->offs[i].uaddr ) ihi = i - 1;
        else if ( uoffset >= fp->idx->offs[i].uaddr ) ilo = i + 1;
        else break;
    }
    int i = ilo-1;
    if (hseek(fp->fp, fp->idx->offs[i].caddr, SEEK_SET) < 0)
    {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }
    fp->block_length = 0;  // indicates current block has not been loaded
    fp->block_address = fp->idx->offs[i].caddr;
    fp->block_offset = 0;
    if ( bgzf_read_block(fp) < 0 ) {
        fp->errcode |= BGZF_ERR_IO;
        return -1;
    }
    if ( uoffset - fp->idx->offs[i].uaddr > 0 )
    {
        fp->block_offset = uoffset - fp->idx->offs[i].uaddr;
        assert( fp->block_offset <= fp->block_length );     // todo: skipped, unindexed, blocks
    }
    fp->uncompressed_address = uoffset;
    return 0;
}

long bgzf_utell(BGZF *fp)
{
    return fp->uncompressed_address;    // currently maintained only when reading
}

/* prototype is in hfile_internal.h */
struct hFILE *bgzf_hfile(struct BGZF *fp) {
    return fp->fp;
}
