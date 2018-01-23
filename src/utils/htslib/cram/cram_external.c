/*
Copyright (c) 2015 Genome Research Ltd.
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

/*! \file
 * External CRAM interface.
 *
 * Internally we're happy to use macros and to grub around in the cram
 * structures.  This isn't very sustainable for an externally usable
 * ABI though, so we have anonymous structs and accessor functions too
 * to permit software such as samtools reheader to manipulate cram
 * containers and blocks in a robust manner.
 */

#include <config.h>

#include "htslib/hfile.h"
#include "cram/cram.h"

/*
 *-----------------------------------------------------------------------------
 * cram_fd
 */
SAM_hdr *cram_fd_get_header(cram_fd *fd) { return fd->header; }
void cram_fd_set_header(cram_fd *fd, SAM_hdr *hdr) { fd->header = hdr; }

int cram_fd_get_version(cram_fd *fd) { return fd->version; }
void cram_fd_set_version(cram_fd *fd, int vers) { fd->version = vers; }

int cram_major_vers(cram_fd *fd) { return CRAM_MAJOR_VERS(fd->version); }
int cram_minor_vers(cram_fd *fd) { return CRAM_MINOR_VERS(fd->version); }

hFILE *cram_fd_get_fp(cram_fd *fd) { return fd->fp; }
void cram_fd_set_fp(cram_fd *fd, hFILE *fp) { fd->fp = fp; }


/*
 *-----------------------------------------------------------------------------
 * cram_container
 */
int32_t cram_container_get_length(cram_container *c) {
    return c->length;
}

void cram_container_set_length(cram_container *c, int32_t length) {
    c->length = length;
}


int32_t cram_container_get_num_blocks(cram_container *c) {
    return c->num_blocks;
}

void cram_container_set_num_blocks(cram_container *c, int32_t num_blocks) {
    c->num_blocks = num_blocks;
}


/* Returns the landmarks[] array and the number of elements
 * in num_landmarks.
 */
int32_t *cram_container_get_landmarks(cram_container *c, int32_t *num_landmarks) {
    *num_landmarks = c->num_landmarks;
    return c->landmark;
}

/* Sets the landmarks[] array (pointer copy, not a memory dup) and
 * num_landmarks value.
 */
void cram_container_set_landmarks(cram_container *c, int32_t num_landmarks,
				  int32_t *landmarks) {
    c->num_landmarks = num_landmarks;
    c->landmark = landmarks;
}


/* Returns true if the container is empty (EOF marker) */
int cram_container_is_empty(cram_fd *fd) {
    return fd->empty_container;
}


/*
 *-----------------------------------------------------------------------------
 * cram_block_compression_hdr
 */

/*
 * Utility function to edit an RG id.
 * This is only possible if there is one single RG value used and it
 * is in the container compression header using HUFFMAN or BETA
 * codec.  In this case it is essentially hard coded and needs no
 * editing of external (or worse, CORE) blocks.
 *
 * Returns 0 on success
 *        -1 on failure
 */
// Or arbitrary set compression header constant?

static int cram_block_compression_hdr_set_DS(cram_block_compression_hdr *ch,
					     int ds, int new_rg) {
    if (!ch || !ch->codecs[ds])
	return -1;

    switch (ch->codecs[ds]->codec) {
    case E_HUFFMAN:
	if (ch->codecs[ds]->huffman.ncodes != 1)
	    return -1;
	ch->codecs[ds]->huffman.codes[0].symbol = new_rg;
	return 0;

    case E_BETA:
	if (ch->codecs[ds]->beta.nbits != 0)
	    return -1;
	ch->codecs[ds]->beta.offset = -new_rg;
	return 0;

    default:
	return -1;
    }

    return 0;
}

int cram_block_compression_hdr_set_rg(cram_block_compression_hdr *ch, int new_rg) {
    return cram_block_compression_hdr_set_DS(ch, DS_RG, new_rg);
}

/*
 * Converts a cram_block_compression_hdr struct used for decoding to
 * one used for encoding.  Maybe this should be a transparent
 * operation applied on-demand.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_block_compression_hdr_decoder2encoder(cram_fd *fd,
					       cram_block_compression_hdr *ch) {
    int i;

    if (!ch)
	return -1;

    for (i = 0; i < DS_END; i++) {
	cram_codec *co = ch->codecs[i];
	if (!co)
	    continue;

	if (-1 == cram_codec_decoder2encoder(fd, co))
	    return -1;
    }

    return 0;
}

/*
 *-----------------------------------------------------------------------------
 * cram_slice
 */
int32_t cram_slice_hdr_get_num_blocks(cram_block_slice_hdr *hdr) {
    return hdr->num_blocks;
}


/*
 *-----------------------------------------------------------------------------
 * cram_block
 */
int32_t cram_block_get_content_id(cram_block *b)  { return b->content_id; }
int32_t cram_block_get_comp_size(cram_block *b)   { return b->comp_size; }
int32_t cram_block_get_uncomp_size(cram_block *b) { return b->uncomp_size; }
int32_t cram_block_get_crc32(cram_block *b)       { return b->crc32; }
void *  cram_block_get_data(cram_block *b)        { return BLOCK_DATA(b); }
int32_t cram_block_get_size(cram_block *b)        { return BLOCK_SIZE(b); }
enum cram_content_type cram_block_get_content_type(cram_block *b) {
    return b->content_type;
}

void cram_block_set_content_id(cram_block *b, int32_t id) { b->content_id = id; }
void cram_block_set_comp_size(cram_block *b, int32_t size) { b->comp_size = size; }
void cram_block_set_uncomp_size(cram_block *b, int32_t size) { b->uncomp_size = size; }
void cram_block_set_crc32(cram_block *b, int32_t crc) { b->crc32 = crc; }
void cram_block_set_data(cram_block *b, void *data) { BLOCK_DATA(b) = data; }
void cram_block_set_size(cram_block *b, int32_t size) { BLOCK_SIZE(b) = size; }

int cram_block_append(cram_block *b, void *data, int size) {
    BLOCK_APPEND(b, data, size);
    return BLOCK_DATA(b) ? 0 : -1; // It'll do for now...
}
void cram_block_update_size(cram_block *b) { BLOCK_UPLEN(b); }

// Offset is known as "size" internally, but it can be confusing.
size_t cram_block_get_offset(cram_block *b) { return BLOCK_SIZE(b); }
void cram_block_set_offset(cram_block *b, size_t offset) { BLOCK_SIZE(b) = offset; }


/*
 * Copies the blocks representing the next num_slice slices from a
 * container from 'in' to 'out'.  It is expected that the file pointer
 * is just after the read of the cram_container and cram compression
 * header.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_copy_slice(cram_fd *in, cram_fd *out, int32_t num_slice) {
    int32_t i, j;

    for (i = 0; i < num_slice; i++) {
	cram_block *blk;
	cram_block_slice_hdr *hdr;

	if (!(blk = cram_read_block(in)))
	    return -1;
	if (!(hdr = cram_decode_slice_header(in, blk))) {
	    cram_free_block(blk);
	    return -1;
	}
	if (cram_write_block(out, blk) != 0) {
	    cram_free_block(blk);
	    return -1;
	}
	cram_free_block(blk);

	int num_blocks = cram_slice_hdr_get_num_blocks(hdr);
	for (j = 0; j < num_blocks; j++) {
	    blk = cram_read_block(in);
	    if (!blk || cram_write_block(out, blk) != 0) {
		if (blk) cram_free_block(blk);
		return -1;
	    }
	    cram_free_block(blk);
	}
	cram_free_slice_header(hdr);
    }

    return 0;
}

/*
 * Renumbers RG numbers in a cram compression header.
 *
 * CRAM stores RG as the Nth number in the header, rather than a
 * string holding the ID: tag.  This is smaller in space, but means
 * "samtools cat" to join files together that contain single but
 * different RG lines needs a way of renumbering them.
 *
 * The file descriptor is expected to be immediately after the
 * cram_container structure (ie before the cram compression header).
 * Due to the nature of the CRAM format, this needs to read and write
 * the blocks itself.  Note that there may be multiple slices within
 * the container, meaning multiple compression headers to manipulate.
 * Changing RG may change the size of the compression header and
 * therefore the length field in the container.  Hence we rewrite all
 * blocks just incase and also emit the adjusted container.
 *
 * The current implementation can only cope with renumbering a single
 * RG (and only then if it is using HUFFMAN or BETA codecs).  In
 * theory it *may* be possible to renumber multiple RGs if they use
 * HUFFMAN to the CORE block or use an external block unshared by any
 * other data series.  So we have an API that can be upgraded to
 * support this, but do not implement it for now.  An example
 * implementation of RG as an EXTERNAL block would be to find that
 * block and rewrite it, returning the number of blocks consumed.
 *
 * Returns 0 on success;
 *        -1 if unable to edit;
 *        -2 on other errors (eg I/O).
 */
int cram_transcode_rg(cram_fd *in, cram_fd *out,
		      cram_container *c,
		      int nrg, int *in_rg, int *out_rg) {
    int new_rg = *out_rg, old_size, new_size;
    cram_block *o_blk, *n_blk;
    cram_block_compression_hdr *ch;

    if (nrg != 1) {
	hts_log_error("CRAM transcode supports only a single RG");
	return -2;
    }

    // Produce a new block holding the updated compression header,
    // with RG transcoded to a new value. (Single only supported.)
    o_blk = cram_read_block(in);
    old_size = cram_block_size(o_blk);
    ch = cram_decode_compression_header(in, o_blk);
    if (cram_block_compression_hdr_set_rg(ch, new_rg) != 0)
	return -1;
    cram_block_compression_hdr_decoder2encoder(in, ch);
    n_blk = cram_encode_compression_header(in, c, ch);
    cram_free_compression_header(ch);

    /*
     * Warning: this has internal knowledge of the cram compression
     * header format.
     *
     * The decoder doesn't set c->tags_used, so the encoder puts a two
     * byte blank segment.  This means n_blk is too short.  We skip
     * through the decoded old block (o_blk) and copy from there.
     */
    char *cp = cram_block_get_data(o_blk);
    char *op = cp;
    char *endp = cp + cram_block_get_uncomp_size(o_blk);
    //fprintf(stderr, "sz = %d\n", (int)(endp-cp));
    int32_t i32;

    cp += safe_itf8_get(cp, endp, &i32);
    cp += i32;
    cp += safe_itf8_get(cp, endp, &i32);
    cp += i32;
    op = cp;
    cp += safe_itf8_get(cp, endp, &i32);
    i32 += (cp-op);

    //fprintf(stderr, "remaining %d bytes\n", i32);
    cram_block_set_size(n_blk, cram_block_get_size(n_blk)-2);
    cram_block_append(n_blk, op, i32);
    cram_block_update_size(n_blk);

    new_size = cram_block_size(n_blk);

    //fprintf(stderr, "size %d -> %d\n", old_size, new_size);

    // Now we've constructedthe updated compression header,
    // amend the container too (it may have changed size).
    int32_t *landmarks, num_landmarks;
    landmarks = cram_container_get_landmarks(c, &num_landmarks);

    if (old_size != new_size) {
	int diff = new_size - old_size, j;

	for (j = 0; j < num_landmarks; j++)
	    landmarks[j] += diff;
	//cram_container_set_landmarks(c, num_landmarks, landmarks);
	cram_container_set_length(c, cram_container_get_length(c) + diff);
    }

    // Finally write it all out; container, compression header,
    // and then all the remaining slice blocks.
    if (cram_write_container(out, c) != 0)
	return -2;

    cram_write_block(out, n_blk);
    cram_free_block(o_blk);
    cram_free_block(n_blk);

    // Container num_blocks can be invalid, due to a bug.
    // Instead we iterate in slice context instead.
    return cram_copy_slice(in, out, num_landmarks);
}


/*!
 * Returns the refs_t structure used by a cram file handle.
 *
 * This may be used in conjunction with option CRAM_OPT_SHARED_REF to
 * share reference memory between multiple file handles.
 *
 * @return
 * Returns NULL if none exists or the file handle is not a CRAM file.
 */
refs_t *cram_get_refs(htsFile *fd) {
    return fd->format.format == cram
        ? fd->fp.cram->refs
        : NULL;
}
