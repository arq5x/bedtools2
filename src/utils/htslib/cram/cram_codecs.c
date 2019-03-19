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

/*
 * FIXME: add checking of cram_external_type to return NULL on unsupported
 * {codec,type} tuples.
 */

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <errno.h>

#include "cram/cram.h"

/*
 * ---------------------------------------------------------------------------
 * Block bit-level I/O functions.
 * All defined static here to promote easy inlining by the compiler.
 */

#if 0
/* Get a single bit, MSB first */
static signed int get_bit_MSB(cram_block *block) {
    unsigned int val;

    if (block->byte > block->alloc)
        return -1;

    val = block->data[block->byte] >> block->bit;
    if (--block->bit == -1) {
        block->bit = 7;
        block->byte++;
        //printf("(%02X)", block->data[block->byte]);
    }

    //printf("-B%d-", val&1);

    return val & 1;
}
#endif

/*
 * Count number of successive 0 and 1 bits
 */
static int get_one_bits_MSB(cram_block *block) {
    int n = 0, b;
    if (block->byte >= block->uncomp_size)
        return -1;
    do {
        b = block->data[block->byte] >> block->bit;
        if (--block->bit == -1) {
            block->bit = 7;
            block->byte++;
            if (block->byte == block->uncomp_size && (b&1))
                return -1;
        }
        n++;
    } while (b&1);

    return n-1;
}

static int get_zero_bits_MSB(cram_block *block) {
    int n = 0, b;
    if (block->byte >= block->uncomp_size)
        return -1;
    do {
        b = block->data[block->byte] >> block->bit;
        if (--block->bit == -1) {
            block->bit = 7;
            block->byte++;
            if (block->byte == block->uncomp_size && !(b&1))
                return -1;
        }
        n++;
    } while (!(b&1));

    return n-1;
}

#if 0
/* Stores a single bit */
static void store_bit_MSB(cram_block *block, unsigned int bit) {
    if (block->byte >= block->alloc) {
        block->alloc = block->alloc ? block->alloc*2 : 1024;
        block->data = realloc(block->data, block->alloc);
    }

    if (bit)
        block->data[block->byte] |= (1 << block->bit);

    if (--block->bit == -1) {
        block->bit = 7;
        block->byte++;
        block->data[block->byte] = 0;
    }
}
#endif

#if 0
/* Rounds to the next whole byte boundary first */
static void store_bytes_MSB(cram_block *block, char *bytes, int len) {
    if (block->bit != 7) {
        block->bit = 7;
        block->byte++;
    }

    while (block->byte + len >= block->alloc) {
        block->alloc = block->alloc ? block->alloc*2 : 1024;
        block->data = realloc(block->data, block->alloc);
    }

    memcpy(&block->data[block->byte], bytes, len);
    block->byte += len;
}
#endif

/* Local optimised copy for inlining */
static inline unsigned int get_bits_MSB(cram_block *block, int nbits) {
    unsigned int val = 0;
    int i;

#if 0
    // Fits within the current byte */
    if (nbits <= block->bit+1) {
        val = (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
        if ((block->bit -= nbits) == -1) {
            block->bit = 7;
            block->byte++;
        }
        return val;
    }

    // partial first byte
    val = block->data[block->byte] & ((1<<(block->bit+1))-1);
    nbits -= block->bit+1;
    block->bit = 7;
    block->byte++;

    // whole middle bytes
    while (nbits >= 8) {
        val = (val << 8) | block->data[block->byte++];
        nbits -= 8;
    }

    val <<= nbits;
    val |= (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
    block->bit -= nbits;
    return val;
#endif

#if 0
    /* Inefficient implementation! */
    //printf("{");
    for (i = 0; i < nbits; i++)
        //val = (val << 1) | get_bit_MSB(block);
        GET_BIT_MSB(block, val);
#endif

#if 1
    /* Combination of 1st two methods */
    if (nbits <= block->bit+1) {
        val = (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
        if ((block->bit -= nbits) == -1) {
            block->bit = 7;
            block->byte++;
        }
        return val;
    }

    switch(nbits) {
//  case 15: GET_BIT_MSB(block, val);
//  case 14: GET_BIT_MSB(block, val);
//  case 13: GET_BIT_MSB(block, val);
//  case 12: GET_BIT_MSB(block, val);
//  case 11: GET_BIT_MSB(block, val);
//  case 10: GET_BIT_MSB(block, val);
//  case  9: GET_BIT_MSB(block, val);
    case  8: GET_BIT_MSB(block, val);
    case  7: GET_BIT_MSB(block, val);
    case  6: GET_BIT_MSB(block, val);
    case  5: GET_BIT_MSB(block, val);
    case  4: GET_BIT_MSB(block, val);
    case  3: GET_BIT_MSB(block, val);
    case  2: GET_BIT_MSB(block, val);
    case  1: GET_BIT_MSB(block, val);
        break;

    default:
        for (i = 0; i < nbits; i++)
            //val = (val << 1) | get_bit_MSB(block);
            GET_BIT_MSB(block, val);
    }
#endif

    //printf("=0x%x}", val);

    return val;
}

/*
 * Can store up to 24-bits worth of data encoded in an integer value
 * Possibly we'd want to have a less optimal store_bits function when dealing
 * with nbits > 24, but for now we assume the codes generated are never
 * that big. (Given this is only possible with 121392 or more
 * characters with exactly the correct frequency distribution we check
 * for it elsewhere.)
 */
static int store_bits_MSB(cram_block *block, unsigned int val, int nbits) {
    //fprintf(stderr, " store_bits: %02x %d\n", val, nbits);

    /*
     * Use slow mode until we tweak the huffman generator to never generate
     * codes longer than 24-bits.
     */
    unsigned int mask;

    if (block->byte+4 >= block->alloc) {
        if (block->byte) {
            block->alloc *= 2;
            block->data = realloc(block->data, block->alloc + 4);
            if (!block->data)
                return -1;
        } else {
            block->alloc = 1024;
            block->data = realloc(block->data, block->alloc + 4);
            if (!block->data)
                return -1;
            block->data[0] = 0; // initialise first byte of buffer
        }
    }

    /* fits in current bit-field */
    if (nbits <= block->bit+1) {
        block->data[block->byte] |= (val << (block->bit+1-nbits));
        if ((block->bit-=nbits) == -1) {
            block->bit = 7;
            block->byte++;
            block->data[block->byte] = 0;
        }
        return 0;
    }

    block->data[block->byte] |= (val >> (nbits -= block->bit+1));
    block->bit = 7;
    block->byte++;
    block->data[block->byte] = 0;

    mask = 1<<(nbits-1);
    do {
        if (val & mask)
            block->data[block->byte] |= (1 << block->bit);
        if (--block->bit == -1) {
            block->bit = 7;
            block->byte++;
            block->data[block->byte] = 0;
        }
        mask >>= 1;
    } while(--nbits);

    return 0;
}

/*
 * Returns the next 'size' bytes from a block, or NULL if insufficient
 * data left.This is just a pointer into the block data and not an
 * allocated object, so do not free the result.
 */
static char *cram_extract_block(cram_block *b, int size) {
    char *cp = (char *)b->data + b->idx;
    b->idx += size;
    if (b->idx > b->uncomp_size)
        return NULL;

    return cp;
}

/*
 * ---------------------------------------------------------------------------
 * EXTERNAL
 */
int cram_external_decode_int(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    int l;
    char *cp;
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    l = safe_itf8_get(cp, (char *)b->data + b->uncomp_size, (int32_t *)out);
    b->idx += l;
    *out_size = 1;

    return l > 0 ? 0 : -1;
}

int cram_external_decode_char(cram_slice *slice, cram_codec *c,
                              cram_block *in, char *out,
                              int *out_size) {
    char *cp;
    cram_block *b;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = cram_extract_block(b, *out_size);
    if (!cp)
        return -1;

    if (out)
        memcpy(out, cp, *out_size);
    return 0;
}

static int cram_external_decode_block(cram_slice *slice, cram_codec *c,
                                      cram_block *in, char *out_,
                                      int *out_size) {
    char *cp;
    cram_block *out = (cram_block *)out_;
    cram_block *b = NULL;

    /* Find the external block */
    b = cram_get_block_by_id(slice, c->external.content_id);
    if (!b)
        return *out_size?-1:0;

    cp = cram_extract_block(b, *out_size);
    if (!cp)
        return -1;

    BLOCK_APPEND(out, cp, *out_size);
    return 0;
}

void cram_external_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

cram_codec *cram_external_decode_init(char *data, int size,
                                      enum cram_external_type option,
                                      int version) {
    cram_codec *c = NULL;
    char *cp = data;

    if (size < 1)
        goto malformed;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_EXTERNAL;
    if (option == E_INT || option == E_LONG)
        c->decode = cram_external_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
        c->decode = cram_external_decode_char;
    else
        c->decode = cram_external_decode_block;
    c->free   = cram_external_decode_free;

    cp += safe_itf8_get(cp, data + size, &c->external.content_id);

    if (cp - data != size)
        goto malformed;

    c->external.type = option;

    return c;

 malformed:
    hts_log_error("Malformed external header stream");
    free(c);
    return NULL;
}

int cram_external_encode_int(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    uint32_t *i32 = (uint32_t *)in;

    itf8_put_blk(c->out, *i32);
    return 0;
}

int cram_external_encode_char(cram_slice *slice, cram_codec *c,
                              char *in, int in_size) {
    BLOCK_APPEND(c->out, in, in_size);
    return 0;
}

void cram_external_encode_free(cram_codec *c) {
    if (!c)
        return;
    free(c);
}

int cram_external_encode_store(cram_codec *c, cram_block *b, char *prefix,
                               int version) {
    char tmp[99], *tp = tmp;
    int len = 0;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    tp += itf8_put(tp, c->e_external.content_id);
    len += itf8_put_blk(b, c->codec);
    len += itf8_put_blk(b, tp-tmp);
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    return len;
}

cram_codec *cram_external_encode_init(cram_stats *st,
                                      enum cram_external_type option,
                                      void *dat,
                                      int version) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec = E_EXTERNAL;
    c->free = cram_external_encode_free;
    if (option == E_INT || option == E_LONG)
        c->encode = cram_external_encode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
        c->encode = cram_external_encode_char;
    else
        abort();
    c->store = cram_external_encode_store;

    c->e_external.content_id = (size_t)dat;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BETA
 */
int cram_beta_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n = *out_size;

    if (c->beta.nbits) {
        if (cram_not_enough_bits(in, c->beta.nbits * n))
            return -1;

        for (i = 0; i < n; i++)
            out_i[i] = get_bits_MSB(in, c->beta.nbits) - c->beta.offset;
    } else {
        for (i = 0; i < n; i++)
            out_i[i] = -c->beta.offset;
    }

    return 0;
}

int cram_beta_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int i, n = *out_size;


    if (c->beta.nbits) {
        if (cram_not_enough_bits(in, c->beta.nbits * n))
            return -1;

        if (out)
            for (i = 0; i < n; i++)
                out[i] = get_bits_MSB(in, c->beta.nbits) - c->beta.offset;
        else
            for (i = 0; i < n; i++)
                get_bits_MSB(in, c->beta.nbits);
    } else {
        if (out)
            for (i = 0; i < n; i++)
                out[i] = -c->beta.offset;
    }

    return 0;
}

void cram_beta_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

cram_codec *cram_beta_decode_init(char *data, int size,
                                  enum cram_external_type option,
                                  int version) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_BETA;
    if (option == E_INT || option == E_LONG)
        c->decode = cram_beta_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
        c->decode = cram_beta_decode_char;
    else {
        hts_log_error("BYTE_ARRAYs not supported by this codec");
        return NULL;
    }
    c->free   = cram_beta_decode_free;

    c->beta.nbits = -1;
    cp += safe_itf8_get(cp, data + size, &c->beta.offset);
    if (cp < data + size) // Ensure test below works
        cp += safe_itf8_get(cp, data + size, &c->beta.nbits);

    if (cp - data != size
        || c->beta.nbits < 0 || c->beta.nbits > 8 * sizeof(int)) {
        hts_log_error("Malformed beta header stream");
        free(c);
        return NULL;
    }

    return c;
}

int cram_beta_encode_store(cram_codec *c, cram_block *b,
                           char *prefix, int version) {
    int len = 0;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    len += itf8_put_blk(b, c->codec);
    len += itf8_put_blk(b, itf8_size(c->e_beta.offset)
                        + itf8_size(c->e_beta.nbits)); // codec length
    len += itf8_put_blk(b, c->e_beta.offset);
    len += itf8_put_blk(b, c->e_beta.nbits);

    return len;
}

int cram_beta_encode_int(cram_slice *slice, cram_codec *c,
                         char *in, int in_size) {
    int *syms = (int *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
        r |= store_bits_MSB(c->out, syms[i] + c->e_beta.offset,
                            c->e_beta.nbits);

    return r;
}

int cram_beta_encode_char(cram_slice *slice, cram_codec *c,
                          char *in, int in_size) {
    unsigned char *syms = (unsigned char *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
        r |= store_bits_MSB(c->out, syms[i] + c->e_beta.offset,
                            c->e_beta.nbits);

    return r;
}

void cram_beta_encode_free(cram_codec *c) {
    if (c) free(c);
}

cram_codec *cram_beta_encode_init(cram_stats *st,
                                  enum cram_external_type option,
                                  void *dat,
                                  int version) {
    cram_codec *c;
    int min_val, max_val, len = 0;
    int64_t range;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec  = E_BETA;
    c->free   = cram_beta_encode_free;
    if (option == E_INT)
        c->encode = cram_beta_encode_int;
    else
        c->encode = cram_beta_encode_char;
    c->store  = cram_beta_encode_store;

    if (dat) {
        min_val = ((int *)dat)[0];
        max_val = ((int *)dat)[1];
    } else {
        min_val = INT_MAX;
        max_val = INT_MIN;
        int i;
        for (i = 0; i < MAX_STAT_VAL; i++) {
            if (!st->freqs[i])
                continue;
            if (min_val > i)
                min_val = i;
            max_val = i;
        }
        if (st->h) {
            khint_t k;

            for (k = kh_begin(st->h); k != kh_end(st->h); k++) {
                if (!kh_exist(st->h, k))
                    continue;

                i = kh_key(st->h, k);
                if (min_val > i)
                    min_val = i;
                if (max_val < i)
                    max_val = i;
            }
        }
    }

    assert(max_val >= min_val);
    c->e_beta.offset = -min_val;
    range = (int64_t) max_val - min_val;
    while (range) {
        len++;
        range >>= 1;
    }
    c->e_beta.nbits = len;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * SUBEXP
 */
int cram_subexp_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int n, count;
    int k = c->subexp.k;

    for (count = 0, n = *out_size; count < n; count++) {
        int i = 0, tail;
        int val;

        /* Get number of 1s */
        //while (get_bit_MSB(in) == 1) i++;
        i = get_one_bits_MSB(in);
        if (i < 0 || cram_not_enough_bits(in, i > 0 ? i + k - 1 : k))
            return -1;
        /*
         * Val is
         * i > 0:  2^(k+i-1) + k+i-1 bits
         * i = 0:  k bits
         */
        if (i) {
            tail = i + k-1;
            val = 0;
            while (tail) {
                //val = val<<1; val |= get_bit_MSB(in);
                GET_BIT_MSB(in, val);
                tail--;
            }
            val += 1 << (i + k-1);
        } else {
            tail = k;
            val = 0;
            while (tail) {
                //val = val<<1; val |= get_bit_MSB(in);
                GET_BIT_MSB(in, val);
                tail--;
            }
        }

        out_i[count] = val - c->subexp.offset;
    }

    return 0;
}

void cram_subexp_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

cram_codec *cram_subexp_decode_init(char *data, int size,
                                    enum cram_external_type option,
                                    int version) {
    cram_codec *c;
    char *cp = data;

    if (option != E_INT) {
        hts_log_error("This codec only supports INT encodings");
        return NULL;
    }

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_SUBEXP;
    c->decode = cram_subexp_decode;
    c->free   = cram_subexp_decode_free;
    c->subexp.k = -1;

    cp += safe_itf8_get(cp, data + size, &c->subexp.offset);
    cp += safe_itf8_get(cp, data + size, &c->subexp.k);

    if (cp - data != size || c->subexp.k < 0) {
        hts_log_error("Malformed subexp header stream");
        free(c);
        return NULL;
    }

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * GAMMA
 */
int cram_gamma_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;

    for (i = 0, n = *out_size; i < n; i++) {
        int nz = 0;
        int val;
        //while (get_bit_MSB(in) == 0) nz++;
        nz = get_zero_bits_MSB(in);
        if (cram_not_enough_bits(in, nz))
            return -1;
        val = 1;
        while (nz > 0) {
            //val <<= 1; val |= get_bit_MSB(in);
            GET_BIT_MSB(in, val);
            nz--;
        }

        out_i[i] = val - c->gamma.offset;
    }

    return 0;
}

void cram_gamma_decode_free(cram_codec *c) {
    if (c)
        free(c);
}

cram_codec *cram_gamma_decode_init(char *data, int size,
                                   enum cram_external_type option,
                                   int version) {
    cram_codec *c = NULL;
    char *cp = data;

    if (option != E_INT) {
        hts_log_error("This codec only supports INT encodings");
        return NULL;
    }

    if (size < 1)
        goto malformed;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_GAMMA;
    c->decode = cram_gamma_decode;
    c->free   = cram_gamma_decode_free;

    cp += safe_itf8_get(cp, data + size, &c->gamma.offset);

    if (cp - data != size)
        goto malformed;

    return c;

 malformed:
    hts_log_error("Malformed gamma header stream");
    free(c);
    return NULL;
}

/*
 * ---------------------------------------------------------------------------
 * HUFFMAN
 */

static int code_sort(const void *vp1, const void *vp2) {
    const cram_huffman_code *c1 = (const cram_huffman_code *)vp1;
    const cram_huffman_code *c2 = (const cram_huffman_code *)vp2;

    if (c1->len != c2->len)
        return c1->len - c2->len;
    else
        return c1->symbol < c2->symbol ? -1 : (c1->symbol > c2->symbol ? 1 : 0);
}

void cram_huffman_decode_free(cram_codec *c) {
    if (!c)
        return;

    if (c->huffman.codes)
        free(c->huffman.codes);
    free(c);
}

int cram_huffman_decode_null(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    return -1;
}

int cram_huffman_decode_char0(cram_slice *slice, cram_codec *c,
                              cram_block *in, char *out, int *out_size) {
    int i, n;

    if (!out)
        return 0;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
        out[i] = c->huffman.codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_char(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    int i, n, ncodes = c->huffman.ncodes;
    const cram_huffman_code * const codes = c->huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
        int idx = 0;
        int val = 0, len = 0, last_len = 0;

        for (;;) {
            int dlen = codes[idx].len - last_len;
            if (cram_not_enough_bits(in, dlen))
                return -1;

            //val <<= dlen;
            //val  |= get_bits_MSB(in, dlen);
            //last_len = (len += dlen);

            last_len = (len += dlen);
            for (; dlen; dlen--) GET_BIT_MSB(in, val);

            idx = val - codes[idx].p;
            if (idx >= ncodes || idx < 0)
                return -1;

            if (codes[idx].code == val && codes[idx].len == len) {
                if (out) out[i] = codes[idx].symbol;
                break;
            }
        }
    }

    return 0;
}

int cram_huffman_decode_int0(cram_slice *slice, cram_codec *c,
                             cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;
    const cram_huffman_code * const codes = c->huffman.codes;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
        out_i[i] = codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_int(cram_slice *slice, cram_codec *c,
                            cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n, ncodes = c->huffman.ncodes;
    const cram_huffman_code * const codes = c->huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
        int idx = 0;
        int val = 0, len = 0, last_len = 0;

        // Now one bit at a time for remaining checks
        for (;;) {
            int dlen = codes[idx].len - last_len;
            if (cram_not_enough_bits(in, dlen))
                return -1;

            //val <<= dlen;
            //val  |= get_bits_MSB(in, dlen);
            //last_len = (len += dlen);

            last_len = (len += dlen);
            for (; dlen; dlen--) GET_BIT_MSB(in, val);

            idx = val - codes[idx].p;
            if (idx >= ncodes || idx < 0)
                return -1;

            if (codes[idx].code == val && codes[idx].len == len) {
                out_i[i] = codes[idx].symbol;
                break;
            }
        }
    }

    return 0;
}

/*
 * Initialises a huffman decoder from an encoding data stream.
 */
cram_codec *cram_huffman_decode_init(char *data, int size,
                                     enum cram_external_type option,
                                     int version) {
    int32_t ncodes = 0, i, j;
    char *cp = data, *data_end = &data[size];
    cram_codec *h;
    cram_huffman_code *codes = NULL;
    int32_t val, last_len, max_len = 0;
    uint32_t max_val; // needs one more bit than val
    const int max_code_bits = sizeof(val) * 8 - 1;
    int l;

    if (option == E_BYTE_ARRAY_BLOCK) {
        hts_log_error("BYTE_ARRAYs not supported by this codec");
        return NULL;
    }

    cp += safe_itf8_get(cp, data_end, &ncodes);
    if (ncodes < 0) {
        hts_log_error("Invalid number of symbols in huffman stream");
        return NULL;
    }
    if (ncodes >= SIZE_MAX / sizeof(*codes)) {
        errno = ENOMEM;
        return NULL;
    }

    h = calloc(1, sizeof(*h));
    if (!h)
        return NULL;

    h->codec  = E_HUFFMAN;
    h->free   = cram_huffman_decode_free;

    h->huffman.ncodes = ncodes;
    if (ncodes) {
        codes = h->huffman.codes = malloc(ncodes * sizeof(*codes));
        if (!codes) {
            free(h);
            return NULL;
        }
    } else {
        codes = h->huffman.codes = NULL;
    }

    /* Read symbols and bit-lengths */
    for (i = 0, l = 1; i < ncodes && l > 0; i++, cp += l) {
        l = safe_itf8_get(cp, data_end, &codes[i].symbol);
    }

    if (l < 1)
        goto malformed;

    cp += safe_itf8_get(cp, data_end, &i);
    if (i != ncodes)
        goto malformed;

    if (ncodes == 0) {
        /* NULL huffman stream.  Ensure it returns an error if
           anything tries to use it. */
        h->decode = cram_huffman_decode_null;
        return h;
    }

    for (i = 0, l = 1; i < ncodes; i++, cp += l) {
        l = safe_itf8_get(cp, data_end, &codes[i].len);
        if (l < 1)
            break;
        if (max_len < codes[i].len)
            max_len = codes[i].len;
    }
    if (l < 1 || cp - data != size || max_len >= ncodes)
        goto malformed;

    /* 31 is max. bits available in val */
    if (max_len > max_code_bits) {
        hts_log_error("Huffman code length (%d) is greater "
                      "than maximum supported (%d)", max_len, max_code_bits);
        free(h);
        free(codes);
        return NULL;
    }

    /* Sort by bit length and then by symbol value */
    qsort(codes, ncodes, sizeof(*codes), code_sort);

    /* Assign canonical codes */
    val = -1, last_len = 0, max_val = 0;
    for (i = 0; i < ncodes; i++) {
        val++;
        if (val > max_val)
            goto malformed;

        if (codes[i].len > last_len) {
            val <<= (codes[i].len - last_len);
            last_len = codes[i].len;
            max_val = (1U << codes[i].len) - 1;
        }
        codes[i].code = val;
    }

    /*
     * Compute the next starting point, offset by the i'th value.
     * For example if codes 10, 11, 12, 13 are 30, 31, 32, 33 then
     * codes[10..13].p = 30 - 10.
     */
    last_len = 0;
    for (i = j = 0; i < ncodes; i++) {
        if (codes[i].len > last_len) {
            j = codes[i].code - i;
            last_len = codes[i].len;
        }
        codes[i].p = j;
    }

    // puts("==HUFF LEN==");
    // for (i = 0; i <= last_len+1; i++) {
    //     printf("len %d=%d prefix %d\n", i, h->huffman.lengths[i], h->huffman.prefix[i]);
    // }
    // puts("===HUFFMAN CODES===");
    // for (i = 0; i < ncodes; i++) {
    //     int j;
    //     printf("%d: %d %d %d ", i, codes[i].symbol, codes[i].len, codes[i].code);
    //     j = codes[i].len;
    //     while (j) {
    //         putchar(codes[i].code & (1 << --j) ? '1' : '0');
    //     }
    //     printf(" %d\n", codes[i].code);
    // }

    if (option == E_BYTE || option == E_BYTE_ARRAY) {
        if (h->huffman.codes[0].len == 0)
            h->decode = cram_huffman_decode_char0;
        else
            h->decode = cram_huffman_decode_char;
    } else if (option == E_BYTE_ARRAY_BLOCK) {
        abort();
    } else {
        if (h->huffman.codes[0].len == 0)
            h->decode = cram_huffman_decode_int0;
        else
            h->decode = cram_huffman_decode_int;
    }

    return (cram_codec *)h;

 malformed:
    hts_log_error("Malformed huffman header stream");
    free(codes);
    free(h);
    return NULL;
}

int cram_huffman_encode_char0(cram_slice *slice, cram_codec *c,
                              char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_char(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    int i, code, len, r = 0;
    unsigned char *syms = (unsigned char *)in;

    while (in_size--) {
        int sym = *syms++;
        if (sym >= -1 && sym < MAX_HUFF) {
            i = c->e_huffman.val2code[sym+1];
            assert(c->e_huffman.codes[i].symbol == sym);
            code = c->e_huffman.codes[i].code;
            len  = c->e_huffman.codes[i].len;
        } else {
            /* Slow - use a lookup table for when sym < MAX_HUFF? */
            for (i = 0; i < c->e_huffman.nvals; i++) {
                if (c->e_huffman.codes[i].symbol == sym)
                    break;
            }
            if (i == c->e_huffman.nvals)
                return -1;

            code = c->e_huffman.codes[i].code;
            len  = c->e_huffman.codes[i].len;
        }

        r |= store_bits_MSB(c->out, code, len);
    }

    return r;
}

int cram_huffman_encode_int0(cram_slice *slice, cram_codec *c,
                             char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_int(cram_slice *slice, cram_codec *c,
                            char *in, int in_size) {
    int i, code, len, r = 0;
    int *syms = (int *)in;

    while (in_size--) {
        int sym = *syms++;

        if (sym >= -1 && sym < MAX_HUFF) {
            i = c->e_huffman.val2code[sym+1];
            assert(c->e_huffman.codes[i].symbol == sym);
            code = c->e_huffman.codes[i].code;
            len  = c->e_huffman.codes[i].len;
        } else {
            /* Slow - use a lookup table for when sym < MAX_HUFFMAN_SYM? */
            for (i = 0; i < c->e_huffman.nvals; i++) {
                if (c->e_huffman.codes[i].symbol == sym)
                    break;
            }
            if (i == c->e_huffman.nvals)
                return -1;

            code = c->e_huffman.codes[i].code;
            len  = c->e_huffman.codes[i].len;
        }

        r |= store_bits_MSB(c->out, code, len);
    }

    return r;
}

void cram_huffman_encode_free(cram_codec *c) {
    if (!c)
        return;

    if (c->e_huffman.codes)
        free(c->e_huffman.codes);
    free(c);
}

/*
 * Encodes a huffman tree.
 * Returns number of bytes written.
 */
int cram_huffman_encode_store(cram_codec *c, cram_block *b, char *prefix,
                              int version) {
    int i, len = 0;
    cram_huffman_code *codes = c->e_huffman.codes;
    /*
     * Up to code length 127 means 2.5e+26 bytes of data required (worst
     * case huffman tree needs symbols with freqs matching the Fibonacci
     * series). So guaranteed 1 byte per code.
     *
     * Symbols themselves could be 5 bytes (eg -1 is 5 bytes in itf8).
     *
     * Therefore 6*ncodes + 5 + 5 + 1 + 5 is max memory
     */
    char *tmp = malloc(6*c->e_huffman.nvals+16);
    char *tp = tmp;

    if (!tmp)
        return -1;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    tp += itf8_put(tp, c->e_huffman.nvals);
    for (i = 0; i < c->e_huffman.nvals; i++) {
        tp += itf8_put(tp, codes[i].symbol);
    }

    tp += itf8_put(tp, c->e_huffman.nvals);
    for (i = 0; i < c->e_huffman.nvals; i++) {
        tp += itf8_put(tp, codes[i].len);
    }

    len += itf8_put_blk(b, c->codec);
    len += itf8_put_blk(b, tp-tmp);
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    free(tmp);

    return len;
}

cram_codec *cram_huffman_encode_init(cram_stats *st,
                                     enum cram_external_type option,
                                     void *dat,
                                     int version) {
    int *vals = NULL, *freqs = NULL, vals_alloc = 0, *lens, code, len;
    int nvals, i, ntot = 0, max_val = 0, min_val = INT_MAX, k;
    cram_codec *c;
    cram_huffman_code *codes;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec = E_HUFFMAN;

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
                free(c);
                return NULL;
            }
        }
        vals[nvals] = i;
        freqs[nvals] = st->freqs[i];
        assert(st->freqs[i] > 0);
        ntot += freqs[nvals];
        if (max_val < i) max_val = i;
        if (min_val > i) min_val = i;
        nvals++;
    }
    if (st->h) {
        khint_t k;

        for (k = kh_begin(st->h); k != kh_end(st->h); k++) {
            if (!kh_exist(st->h, k))
                continue;
            if (nvals >= vals_alloc) {
                vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
                vals  = realloc(vals,  vals_alloc * sizeof(int));
                freqs = realloc(freqs, vals_alloc * sizeof(int));
                if (!vals || !freqs)
                    return NULL;
            }
            vals[nvals]= kh_key(st->h, k);
            freqs[nvals] = kh_val(st->h, k);
            assert(freqs[nvals] > 0);
            ntot += freqs[nvals];
            if (max_val < i) max_val = i;
            if (min_val > i) min_val = i;
            nvals++;
        }
    }

    assert(nvals > 0);

    freqs = realloc(freqs, 2*nvals*sizeof(*freqs));
    lens = calloc(2*nvals, sizeof(*lens));
    if (!lens || !freqs)
        return NULL;

    /* Inefficient, use pointers to form chain so we can insert and maintain
     * a sorted list? This is currently O(nvals^2) complexity.
     */
    for (;;) {
        int low1 = INT_MAX, low2 = INT_MAX;
        int ind1 = 0, ind2 = 0;
        for (i = 0; i < nvals; i++) {
            if (freqs[i] < 0)
                continue;
            if (low1 > freqs[i])
                low2 = low1, ind2 = ind1, low1 = freqs[i], ind1 = i;
            else if (low2 > freqs[i])
                low2 = freqs[i], ind2 = i;
        }
        if (low2 == INT_MAX)
            break;

        freqs[nvals] = low1 + low2;
        lens[ind1] = nvals;
        lens[ind2] = nvals;
        freqs[ind1] *= -1;
        freqs[ind2] *= -1;
        nvals++;
    }
    nvals = nvals/2+1;

    /* Assign lengths */
    for (i = 0; i < nvals; i++) {
        int code_len = 0;
        for (k = lens[i]; k; k = lens[k])
            code_len++;
        lens[i] = code_len;
        freqs[i] *= -1;
        //fprintf(stderr, "%d / %d => %d\n", vals[i], freqs[i], lens[i]);
    }


    /* Sort, need in a struct */
    if (!(codes = malloc(nvals * sizeof(*codes))))
        return NULL;
    for (i = 0; i < nvals; i++) {
        codes[i].symbol = vals[i];
        codes[i].len = lens[i];
    }
    qsort(codes, nvals, sizeof(*codes), code_sort);

    /*
     * Generate canonical codes from lengths.
     * Sort by length.
     * Start with 0.
     * Every new code of same length is +1.
     * Every new code of new length is +1 then <<1 per extra length.
     *
     * /\
     * a/\
     * /\/\
     * bcd/\
     *    ef
     *
     * a 1  0
     * b 3  4 (0+1)<<2
     * c 3  5
     * d 3  6
     * e 4  14  (6+1)<<1
     * f 5  15
     */
    code = 0; len = codes[0].len;
    for (i = 0; i < nvals; i++) {
        while (len != codes[i].len) {
            code<<=1;
            len++;
        }
        codes[i].code = code++;

        if (codes[i].symbol >= -1 && codes[i].symbol < MAX_HUFF)
            c->e_huffman.val2code[codes[i].symbol+1] = i;

        //fprintf(stderr, "sym %d, code %d, len %d\n",
        //      codes[i].symbol, codes[i].code, codes[i].len);
    }

    free(lens);
    free(vals);
    free(freqs);

    c->e_huffman.codes = codes;
    c->e_huffman.nvals = nvals;

    c->free = cram_huffman_encode_free;
    if (option == E_BYTE || option == E_BYTE_ARRAY) {
        if (c->e_huffman.codes[0].len == 0)
            c->encode = cram_huffman_encode_char0;
        else
            c->encode = cram_huffman_encode_char;
    } else {
        if (c->e_huffman.codes[0].len == 0)
            c->encode = cram_huffman_encode_int0;
        else
            c->encode = cram_huffman_encode_int;
    }
    c->store = cram_huffman_encode_store;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_LEN
 */
int cram_byte_array_len_decode(cram_slice *slice, cram_codec *c,
                               cram_block *in, char *out,
                               int *out_size) {
    /* Fetch length */
    int32_t len = 0, one = 1;
    int r;

    r = c->byte_array_len.len_codec->decode(slice, c->byte_array_len.len_codec,
                                            in, (char *)&len, &one);
    //printf("ByteArray Len=%d\n", len);

    if (!r && c->byte_array_len.val_codec && len >= 0) {
        r = c->byte_array_len.val_codec->decode(slice,
                                                c->byte_array_len.val_codec,
                                                in, out, &len);
    } else {
        return -1;
    }

    *out_size = len;

    return r;
}

void cram_byte_array_len_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->byte_array_len.len_codec)
        c->byte_array_len.len_codec->free(c->byte_array_len.len_codec);

    if (c->byte_array_len.val_codec)
        c->byte_array_len.val_codec->free(c->byte_array_len.val_codec);

    free(c);
}

cram_codec *cram_byte_array_len_decode_init(char *data, int size,
                                            enum cram_external_type option,
                                            int version) {
    cram_codec *c;
    char *cp   = data;
    char *endp = data + size;
    int32_t encoding = 0;
    int32_t sub_size = -1;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_BYTE_ARRAY_LEN;
    c->decode = cram_byte_array_len_decode;
    c->free   = cram_byte_array_len_decode_free;

    cp += safe_itf8_get(cp, endp, &encoding);
    cp += safe_itf8_get(cp, endp, &sub_size);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->byte_array_len.len_codec = cram_decoder_init(encoding, cp, sub_size,
                                                    E_INT, version);
    if (c->byte_array_len.len_codec == NULL)
        goto no_codec;
    cp += sub_size;

    sub_size = -1;
    cp += safe_itf8_get(cp, endp, &encoding);
    cp += safe_itf8_get(cp, endp, &sub_size);
    if (sub_size < 0 || endp - cp < sub_size)
        goto malformed;
    c->byte_array_len.val_codec = cram_decoder_init(encoding, cp, sub_size,
                                                    option, version);
    if (c->byte_array_len.val_codec == NULL)
        goto no_codec;
    cp += sub_size;

    if (cp - data != size)
        goto malformed;

    return c;

 malformed:
    hts_log_error("Malformed byte_array_len header stream");
 no_codec:
    free(c);
    return NULL;
}

int cram_byte_array_len_encode(cram_slice *slice, cram_codec *c,
                               char *in, int in_size) {
    int32_t i32 = in_size;
    int r = 0;

    r |= c->e_byte_array_len.len_codec->encode(slice,
                                               c->e_byte_array_len.len_codec,
                                               (char *)&i32, 1);
    r |= c->e_byte_array_len.val_codec->encode(slice,
                                               c->e_byte_array_len.val_codec,
                                               in, in_size);
    return r;
}

void cram_byte_array_len_encode_free(cram_codec *c) {
    if (!c)
        return;

    if (c->e_byte_array_len.len_codec)
        c->e_byte_array_len.len_codec->free(c->e_byte_array_len.len_codec);

    if (c->e_byte_array_len.val_codec)
        c->e_byte_array_len.val_codec->free(c->e_byte_array_len.val_codec);

    free(c);
}

int cram_byte_array_len_encode_store(cram_codec *c, cram_block *b,
                                     char *prefix, int version) {
    int len = 0, len2, len3;
    cram_codec *tc;
    cram_block *b_len, *b_val;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    tc = c->e_byte_array_len.len_codec;
    b_len = cram_new_block(0, 0);
    len2 = tc->store(tc, b_len, NULL, version);

    tc = c->e_byte_array_len.val_codec;
    b_val = cram_new_block(0, 0);
    len3 = tc->store(tc, b_val, NULL, version);

    len += itf8_put_blk(b, c->codec);
    len += itf8_put_blk(b, len2+len3);
    BLOCK_APPEND(b, BLOCK_DATA(b_len), BLOCK_SIZE(b_len));
    BLOCK_APPEND(b, BLOCK_DATA(b_val), BLOCK_SIZE(b_val));

    cram_free_block(b_len);
    cram_free_block(b_val);

    return len + len2 + len3;
}

cram_codec *cram_byte_array_len_encode_init(cram_stats *st,
                                            enum cram_external_type option,
                                            void *dat,
                                            int version) {
    cram_codec *c;
    cram_byte_array_len_encoder *e = (cram_byte_array_len_encoder *)dat;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec = E_BYTE_ARRAY_LEN;
    c->free = cram_byte_array_len_encode_free;
    c->encode = cram_byte_array_len_encode;
    c->store = cram_byte_array_len_encode_store;

    c->e_byte_array_len.len_codec = cram_encoder_init(e->len_encoding,
                                                      st, E_INT,
                                                      e->len_dat,
                                                      version);
    c->e_byte_array_len.val_codec = cram_encoder_init(e->val_encoding,
                                                      NULL, E_BYTE_ARRAY,
                                                      e->val_dat,
                                                      version);

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_STOP
 */
static int cram_byte_array_stop_decode_char(cram_slice *slice, cram_codec *c,
                                            cram_block *in, char *out,
                                            int *out_size) {
    char *cp, ch;
    cram_block *b = NULL;

    b = cram_get_block_by_id(slice, c->byte_array_stop.content_id);
    if (!b)
        return *out_size?-1:0;

    if (b->idx >= b->uncomp_size)
        return -1;

    cp = (char *)b->data + b->idx;
    if (out) {
        while ((ch = *cp) != (char)c->byte_array_stop.stop) {
            if (cp - (char *)b->data >= b->uncomp_size)
                return -1;
            *out++ = ch;
            cp++;
        }
    } else {
        // Consume input, but produce no output
        while ((ch = *cp) != (char)c->byte_array_stop.stop) {
            if (cp - (char *)b->data >= b->uncomp_size)
                return -1;
            cp++;
        }
    }

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

int cram_byte_array_stop_decode_block(cram_slice *slice, cram_codec *c,
                                      cram_block *in, char *out_,
                                      int *out_size) {
    cram_block *b;
    cram_block *out = (cram_block *)out_;
    char *cp, *out_cp, *cp_end;
    char stop;

    b = cram_get_block_by_id(slice, c->byte_array_stop.content_id);
    if (!b)
        return *out_size?-1:0;

    if (b->idx >= b->uncomp_size)
        return -1;
    cp = (char *)b->data + b->idx;
    cp_end = (char *)b->data + b->uncomp_size;
    out_cp = (char *)BLOCK_END(out);

    stop = c->byte_array_stop.stop;
    if (cp_end - cp < out->alloc - out->byte) {
        while (cp != cp_end && *cp != stop)
            *out_cp++ = *cp++;
        BLOCK_SIZE(out) = out_cp - (char *)BLOCK_DATA(out);
    } else {
        char *cp_start;
        for (cp_start = cp; cp != cp_end && *cp != stop; cp++)
            ;
        BLOCK_APPEND(out, cp_start, cp - cp_start);
        BLOCK_GROW(out, cp - cp_start);
    }

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

void cram_byte_array_stop_decode_free(cram_codec *c) {
    if (!c) return;

    free(c);
}

cram_codec *cram_byte_array_stop_decode_init(char *data, int size,
                                             enum cram_external_type option,
                                             int version) {
    cram_codec *c = NULL;
    unsigned char *cp = (unsigned char *)data;

    if (size < (CRAM_MAJOR_VERS(version) == 1 ? 5 : 2))
        goto malformed;

    if (!(c = malloc(sizeof(*c))))
        return NULL;

    c->codec  = E_BYTE_ARRAY_STOP;
    switch (option) {
    case E_BYTE_ARRAY_BLOCK:
        c->decode = cram_byte_array_stop_decode_block;
        break;
    case E_BYTE_ARRAY:
        c->decode = cram_byte_array_stop_decode_char;
        break;
    default:
        hts_log_error("The byte_array_stop codec only supports BYTE_ARRAYs");
        free(c);
        return NULL;
    }
    c->free   = cram_byte_array_stop_decode_free;

    c->byte_array_stop.stop = *cp++;
    if (CRAM_MAJOR_VERS(version) == 1) {
        c->byte_array_stop.content_id = cp[0] + (cp[1]<<8) + (cp[2]<<16)
            + (cp[3]<<24);
        cp += 4;
    } else {
        cp += safe_itf8_get((char *) cp, data + size,
                            &c->byte_array_stop.content_id);
    }

    if ((char *)cp - data != size)
        goto malformed;

    return c;

 malformed:
    hts_log_error("Malformed byte_array_stop header stream");
    free(c);
    return NULL;
}

int cram_byte_array_stop_encode(cram_slice *slice, cram_codec *c,
                                char *in, int in_size) {
    BLOCK_APPEND(c->out, in, in_size);
    BLOCK_APPEND_CHAR(c->out, c->e_byte_array_stop.stop);
    return 0;
}

void cram_byte_array_stop_encode_free(cram_codec *c) {
    if (!c)
        return;
    free(c);
}

int cram_byte_array_stop_encode_store(cram_codec *c, cram_block *b,
                                      char *prefix, int version) {
    int len = 0;
    char buf[20], *cp = buf;

    if (prefix) {
        size_t l = strlen(prefix);
        BLOCK_APPEND(b, prefix, l);
        len += l;
    }

    cp += itf8_put(cp, c->codec);

    if (CRAM_MAJOR_VERS(version) == 1) {
        cp += itf8_put(cp, 5);
        *cp++ = c->e_byte_array_stop.stop;
        *cp++ = (c->e_byte_array_stop.content_id >>  0) & 0xff;
        *cp++ = (c->e_byte_array_stop.content_id >>  8) & 0xff;
        *cp++ = (c->e_byte_array_stop.content_id >> 16) & 0xff;
        *cp++ = (c->e_byte_array_stop.content_id >> 24) & 0xff;
    } else {
        cp += itf8_put(cp, 1 + itf8_size(c->e_byte_array_stop.content_id));
        *cp++ = c->e_byte_array_stop.stop;
        cp += itf8_put(cp, c->e_byte_array_stop.content_id);
    }

    BLOCK_APPEND(b, buf, cp-buf);
    len += cp-buf;

    return len;
}

cram_codec *cram_byte_array_stop_encode_init(cram_stats *st,
                                             enum cram_external_type option,
                                             void *dat,
                                             int version) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
        return NULL;
    c->codec = E_BYTE_ARRAY_STOP;
    c->free = cram_byte_array_stop_encode_free;
    c->encode = cram_byte_array_stop_encode;
    c->store = cram_byte_array_stop_encode_store;

    c->e_byte_array_stop.stop = ((int *)dat)[0];
    c->e_byte_array_stop.content_id = ((int *)dat)[1];

    return c;
}

/*
 * ---------------------------------------------------------------------------
 */

const char *cram_encoding2str(enum cram_encoding t) {
    switch (t) {
    case E_NULL:            return "NULL";
    case E_EXTERNAL:        return "EXTERNAL";
    case E_GOLOMB:          return "GOLOMB";
    case E_HUFFMAN:         return "HUFFMAN";
    case E_BYTE_ARRAY_LEN:  return "BYTE_ARRAY_LEN";
    case E_BYTE_ARRAY_STOP: return "BYTE_ARRAY_STOP";
    case E_BETA:            return "BETA";
    case E_SUBEXP:          return "SUBEXP";
    case E_GOLOMB_RICE:     return "GOLOMB_RICE";
    case E_GAMMA:           return "GAMMA";
    case E_NUM_CODECS:
    default:                return "?";
    }
}

static cram_codec *(*decode_init[])(char *data,
                                    int size,
                                    enum cram_external_type option,
                                    int version) = {
    NULL,
    cram_external_decode_init,
    NULL,
    cram_huffman_decode_init,
    cram_byte_array_len_decode_init,
    cram_byte_array_stop_decode_init,
    cram_beta_decode_init,
    cram_subexp_decode_init,
    NULL,
    cram_gamma_decode_init,
};

cram_codec *cram_decoder_init(enum cram_encoding codec,
                              char *data, int size,
                              enum cram_external_type option,
                              int version) {
    if (codec >= E_NULL && codec < E_NUM_CODECS && decode_init[codec]) {
        return decode_init[codec](data, size, option, version);
    } else {
        hts_log_error("Unimplemented codec of type %s", cram_encoding2str(codec));
        return NULL;
    }
}

static cram_codec *(*encode_init[])(cram_stats *stx,
                                    enum cram_external_type option,
                                    void *opt,
                                    int version) = {
    NULL,
    cram_external_encode_init,
    NULL,
    cram_huffman_encode_init,
    cram_byte_array_len_encode_init,
    cram_byte_array_stop_encode_init,
    cram_beta_encode_init,
    NULL, //cram_subexp_encode_init,
    NULL,
    NULL, //cram_gamma_encode_init,
};

cram_codec *cram_encoder_init(enum cram_encoding codec,
                              cram_stats *st,
                              enum cram_external_type option,
                              void *dat,
                              int version) {
    if (st && !st->nvals)
        return NULL;

    if (encode_init[codec]) {
        cram_codec *r;
        if ((r = encode_init[codec](st, option, dat, version)))
            r->out = NULL;
        return r;
    } else {
        hts_log_error("Unimplemented codec of type %s", cram_encoding2str(codec));
        abort();
    }
}

/*
 * Returns the content_id used by this codec, also in id2 if byte_array_len.
 * Returns -1 for the CORE block and -2 for unneeded.
 * id2 is only filled out for BYTE_ARRAY_LEN which uses 2 codecs.
 */
int cram_codec_to_id(cram_codec *c, int *id2) {
    int bnum1, bnum2 = -2;

    switch (c->codec) {
    case E_HUFFMAN:
        bnum1 = c->huffman.ncodes == 1 ? -2 : -1;
        break;
    case E_GOLOMB:
    case E_BETA:
    case E_SUBEXP:
    case E_GOLOMB_RICE:
    case E_GAMMA:
        bnum1 = -1;
        break;
    case E_EXTERNAL:
        bnum1 = c->external.content_id;
        break;
    case E_BYTE_ARRAY_LEN:
        bnum1 = cram_codec_to_id(c->byte_array_len.len_codec, NULL);
        bnum2 = cram_codec_to_id(c->byte_array_len.val_codec, NULL);
        break;
    case E_BYTE_ARRAY_STOP:
        bnum1 = c->byte_array_stop.content_id;
        break;
    case E_NULL:
        bnum1 = -2;
        break;
    default:
        hts_log_error("Unknown codec type %d", c->codec);
        bnum1 = -1;
    }

    if (id2)
        *id2 = bnum2;
    return bnum1;
}


/*
 * cram_codec structures are specialised for decoding or encoding.
 * Unfortunately this makes turning a decoder into an encoder (such as
 * when transcoding files) problematic.
 *
 * This function converts a cram decoder codec into an encoder version
 * in-place (ie it modifiers the codec itself).
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int cram_codec_decoder2encoder(cram_fd *fd, cram_codec *c) {
    int j;

    switch (c->codec) {
    case E_EXTERNAL:
        // shares struct with decode
        c->free = cram_external_encode_free;
        c->store = cram_external_encode_store;
        if (c->decode == cram_external_decode_int)
            c->encode = cram_external_encode_int;
        else if (c->decode == cram_external_decode_char)
            c->encode = cram_external_encode_char;
        else
            return -1;
        break;

    case E_HUFFMAN: {
        // New structure, so switch.
        // FIXME: we huffman and e_huffman structs amended, we could
        // unify this.
        cram_codec *t = malloc(sizeof(*t));
        t->codec = E_HUFFMAN;
        t->free = cram_huffman_encode_free;
        t->store = cram_huffman_encode_store;
        t->e_huffman.codes = c->huffman.codes;
        t->e_huffman.nvals = c->huffman.ncodes;
        for (j = 0; j < t->e_huffman.nvals; j++) {
            int32_t sym = t->e_huffman.codes[j].symbol;
            if (sym >= -1 && sym < MAX_HUFF)
                t->e_huffman.val2code[sym+1] = j;
        }

        if (c->decode == cram_huffman_decode_char0)
            t->encode = cram_huffman_encode_char0;
        else if (c->decode == cram_huffman_decode_char)
            t->encode = cram_huffman_encode_char;
        else if (c->decode == cram_huffman_decode_int0)
            t->encode = cram_huffman_encode_int0;
        else if (c->decode == cram_huffman_decode_int)
            t->encode = cram_huffman_encode_int;
        else {
            free(t);
            return -1;
        }
        *c = *t;
        free(t);
        break;
    }

    case E_BETA:
        // shares struct with decode
        c->free = cram_beta_encode_free;
        c->store = cram_beta_encode_store;
        if (c->decode == cram_beta_decode_int)
            c->encode = cram_beta_encode_int;
        else if (c->decode == cram_beta_decode_char)
            c->encode = cram_beta_encode_char;
        else
            return -1;
        break;

    case E_BYTE_ARRAY_LEN: {
        cram_codec *t = malloc(sizeof(*t));
        t->codec = E_BYTE_ARRAY_LEN;
        t->free   = cram_byte_array_len_encode_free;
        t->store  = cram_byte_array_len_encode_store;
        t->encode = cram_byte_array_len_encode;
        t->e_byte_array_len.len_codec = c->byte_array_len.len_codec;
        t->e_byte_array_len.val_codec = c->byte_array_len.val_codec;
        if (cram_codec_decoder2encoder(fd, t->e_byte_array_len.len_codec) == -1 ||
            cram_codec_decoder2encoder(fd, t->e_byte_array_len.val_codec) == -1) {
            t->free(t);
            return -1;
        }

        // {len,val}_{encoding,dat} are undefined, but unused.
        // Leaving them unset here means we can test that assertion.
        *c = *t;
        free(t);
        break;
    }

    case E_BYTE_ARRAY_STOP:
        // shares struct with decode
        c->free   = cram_byte_array_stop_encode_free;
        c->store  = cram_byte_array_stop_encode_store;
        c->encode = cram_byte_array_stop_encode;
        break;

    default:
        return -1;
    }

    return 0;
}
