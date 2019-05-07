#ifndef LZMA_STUB_H
#define LZMA_STUB_H

/* Some platforms, notably macOS, ship a usable liblzma shared library but
   do not ship any LZMA header files. The <lzma.h> and <lzma/{*}.h> header
   files that come with the library contain the following statement:

     *
     * Author: Lasse Collin
     *
     * This file has been put into the public domain.
     * You can do whatever you want with this file.
     *

   Accordingly the following declarations have been copied and distilled
   from <lzma/base.h> and <lzma/container.h> (primarily) and are sufficient
   to compile cram/cram_io.c in the absence of proper LZMA headers.

   This file, lzma_stub.h, remains in the public domain.  */

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { LZMA_OK = 0, LZMA_STREAM_END = 1 } lzma_ret;
typedef enum { LZMA_RUN = 0, LZMA_FINISH = 3 } lzma_action;
typedef enum { LZMA_CHECK_CRC32 = 1 } lzma_check;
typedef enum { LZMA_RESERVED_ENUM = 0 } lzma_reserved_enum;

struct lzma_allocator;
struct lzma_internal;

typedef struct {
    const uint8_t *next_in;
    size_t avail_in;
    uint64_t total_in;

    uint8_t *next_out;
    size_t avail_out;
    uint64_t total_out;

    const struct lzma_allocator *allocator;
    struct lzma_internal *internal;

    void *reserved_ptr1;
    void *reserved_ptr2;
    void *reserved_ptr3;
    void *reserved_ptr4;
    uint64_t reserved_int1;
    uint64_t reserved_int2;
    size_t reserved_int3;
    size_t reserved_int4;
    lzma_reserved_enum reserved_enum1;
    lzma_reserved_enum reserved_enum2;
} lzma_stream;

#define LZMA_STREAM_INIT \
    { NULL, 0, 0, NULL, 0, 0, NULL, NULL, \
    NULL, NULL, NULL, NULL, 0, 0, 0, 0, \
    LZMA_RESERVED_ENUM, LZMA_RESERVED_ENUM }

extern size_t lzma_stream_buffer_bound(size_t uncompressed_size);

extern lzma_ret lzma_easy_buffer_encode(
        uint32_t preset, lzma_check check,
        const struct lzma_allocator *allocator,
        const uint8_t *in, size_t in_size,
        uint8_t *out, size_t *out_pos, size_t out_size);

extern lzma_ret lzma_stream_decoder(
        lzma_stream *strm, uint64_t memlimit, uint32_t flags);

extern uint64_t lzma_easy_decoder_memusage(uint32_t preset);

extern lzma_ret lzma_code(lzma_stream *strm, lzma_action action);

extern void lzma_end(lzma_stream *strm);

#ifdef __cplusplus
}
#endif

#endif
