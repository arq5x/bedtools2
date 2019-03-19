/// @file htslib/hfile.h
/// Buffered low-level input/output streams.
/*
    Copyright (C) 2013-2016 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

#ifndef HTSLIB_HFILE_H
#define HTSLIB_HFILE_H

#include <string.h>

#include <sys/types.h>

#include "hts_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct hFILE_backend;
/// Low-level input/output stream handle
/** The fields of this structure are declared here solely for the benefit
of the hFILE-related inline functions.  They may change in future releases.
User code should not use them directly; you should imagine that hFILE is an
opaque incomplete type.
*/
typedef struct hFILE {
    // @cond internal
    char *buffer, *begin, *end, *limit;
    const struct hFILE_backend *backend;
    off_t offset;
    unsigned at_eof:1, mobile:1, readonly:1;
    int has_errno;
    // @endcond
} hFILE;

/// Defines the operations that used by an IO file
typedef struct  hFILE_callback_ops {
	void* cb_data;
	ssize_t (*read)(void* cb_data, void* buf, size_t sz);
	ssize_t (*write)(void* cb_data, const void* buf, size_t sz);
	off_t (*seek)(void* cb_data, off_t ofs, int whence);
	int (*flush)(void* cb_data);
	int (*close)(void* cb_data);
} hFILE_callback_ops;

/// Open the named file or URL as a stream
/** @return An hFILE pointer, or `NULL` (with _errno_ set) if an error occurred.

The usual `fopen(3)` _mode_ letters are supported: one of
`r` (read), `w` (write), `a` (append), optionally followed by any of
`+` (update), `e` (close on `exec(2)`), `x` (create exclusively),
`:` (indicates scheme-specific variable arguments follow).
*/
hFILE *hopen(const char *filename, const char *mode, ...) HTS_RESULT_USED;

/// Wrap a group of operation callbacks into a HFile stream
/** @return An hFILE pointer, or NULL if error occurred
 */
hFILE *hopen_callback(hFILE_callback_ops ops, const char* mode) HTS_RESULT_USED;

/// Associate a stream with an existing open file descriptor
/** @return An hFILE pointer, or `NULL` (with _errno_ set) if an error occurred.

Note that the file must be opened in binary mode, or else
there will be problems on platforms that make a difference
between text and binary mode.

For socket descriptors (on Windows), _mode_ should contain `s`.
*/
hFILE *hdopen(int fd, const char *mode) HTS_RESULT_USED;

/// Report whether the file name or URL denotes remote storage
/** @return  0 if local, 1 if remote.

"Remote" means involving e.g. explicit network access, with the implication
that callers may wish to cache such files' contents locally.
*/
int hisremote(const char *filename) HTS_RESULT_USED;

/// Flush (for output streams) and close the stream
/** @return  0 if successful, or `EOF` (with _errno_ set) if an error occurred.
*/
int hclose(hFILE *fp) HTS_RESULT_USED;

/// Close the stream, without flushing or propagating errors
/** For use while cleaning up after an error only.  Preserves _errno_.
*/
void hclose_abruptly(hFILE *fp);

/// Return the stream's error indicator
/** @return  Non-zero (in fact, an _errno_ value) if an error has occurred.

This would be called `herror()` and return true/false to parallel `ferror(3)`,
but a networking-related `herror(3)` function already exists.
*/
static inline int herrno(hFILE *fp)
{
    return fp->has_errno;
}

/// Clear the stream's error indicator
static inline void hclearerr(hFILE *fp)
{
    fp->has_errno = 0;
}

/// Reposition the read/write stream offset
/** @return  The resulting offset within the stream (as per `lseek(2)`),
    or negative if an error occurred.
*/
off_t hseek(hFILE *fp, off_t offset, int whence) HTS_RESULT_USED;

/// Report the current stream offset
/** @return  The offset within the stream, starting from zero.
*/
static inline off_t htell(hFILE *fp)
{
    return fp->offset + (fp->begin - fp->buffer);
}

/// Read one character from the stream
/** @return  The character read, or `EOF` on end-of-file or error.
*/
static inline int hgetc(hFILE *fp)
{
    extern int hgetc2(hFILE *);
    return (fp->end > fp->begin)? (unsigned char) *(fp->begin++) : hgetc2(fp);
}

/// Read from the stream until the delimiter, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer
    @param delim   The delimiter (interpreted as an `unsigned char`)
    @param fp      The file stream
    @return  The number of bytes read, or negative on error.
    @since   1.4

Bytes will be read into the buffer up to and including a delimiter, until
EOF is reached, or _size-1_ bytes have been written, whichever comes first.
The string will then be terminated with a NUL byte (`\0`).
*/
ssize_t hgetdelim(char *buffer, size_t size, int delim, hFILE *fp)
    HTS_RESULT_USED;

/// Read a line from the stream, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer
    @param fp      The file stream
    @return  The number of bytes read, or negative on error.
    @since   1.4

Specialization of hgetdelim() for a `\n` delimiter.
*/
static inline ssize_t HTS_RESULT_USED
hgetln(char *buffer, size_t size, hFILE *fp)
{
    return hgetdelim(buffer, size, '\n', fp);
}

/// Read a line from the stream, up to a maximum length
/** @param buffer  The buffer into which bytes will be written
    @param size    The size of the buffer (must be > 1 to be useful)
    @param fp      The file stream
    @return  _buffer_ on success, or `NULL` if an error occurred.
    @since   1.4

This function can be used as a replacement for `fgets(3)`, or together with
kstring's `kgetline()` to read arbitrarily-long lines into a _kstring_t_.
*/
char *hgets(char *buffer, int size, hFILE *fp) HTS_RESULT_USED;

/// Peek at characters to be read without removing them from buffers
/** @param fp      The file stream
    @param buffer  The buffer to which the peeked bytes will be written
    @param nbytes  The number of bytes to peek at; limited by the size of the
                   internal buffer, which could be as small as 4K.
    @return  The number of bytes peeked, which may be less than _nbytes_
             if EOF is encountered; or negative, if there was an I/O error.

The characters peeked at remain in the stream's internal buffer, and will be
returned by later hread() etc calls.
*/
ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes) HTS_RESULT_USED;

/// Read a block of characters from the file
/** @return  The number of bytes read, or negative if an error occurred.

The full _nbytes_ requested will be returned, except as limited by EOF
or I/O errors.
*/
static inline ssize_t HTS_RESULT_USED
hread(hFILE *fp, void *buffer, size_t nbytes)
{
    extern ssize_t hread2(hFILE *, void *, size_t, size_t);

    size_t n = fp->end - fp->begin;
    if (n > nbytes) n = nbytes;
    memcpy(buffer, fp->begin, n);
    fp->begin += n;
    return (n == nbytes || !fp->mobile)? (ssize_t) n : hread2(fp, buffer, nbytes, n);
}

/// Write a character to the stream
/** @return  The character written, or `EOF` if an error occurred.
*/
static inline int hputc(int c, hFILE *fp)
{
    extern int hputc2(int, hFILE *);
    if (fp->begin < fp->limit) *(fp->begin++) = c;
    else c = hputc2(c, fp);
    return c;
}

/// Write a string to the stream
/** @return  0 if successful, or `EOF` if an error occurred.
*/
static inline int hputs(const char *text, hFILE *fp)
{
    extern int hputs2(const char *, size_t, size_t, hFILE *);

    size_t nbytes = strlen(text), n = fp->limit - fp->begin;
    if (n > nbytes) n = nbytes;
    memcpy(fp->begin, text, n);
    fp->begin += n;
    return (n == nbytes)? 0 : hputs2(text, nbytes, n, fp);
}

/// Write a block of characters to the file
/** @return  Either _nbytes_, or negative if an error occurred.

In the absence of I/O errors, the full _nbytes_ will be written.
*/
static inline ssize_t HTS_RESULT_USED
hwrite(hFILE *fp, const void *buffer, size_t nbytes)
{
    extern ssize_t hwrite2(hFILE *, const void *, size_t, size_t);
    extern int hfile_set_blksize(hFILE *fp, size_t bufsiz);

    if(!fp->mobile){
        if (fp->limit - fp->begin < (ssize_t)nbytes){
            hfile_set_blksize(fp, fp->limit - fp->buffer + nbytes);
            fp->end = fp->limit;
        }
    }

    size_t n = fp->limit - fp->begin;
    if (nbytes >= n && fp->begin == fp->buffer) {
        // Go straight to hwrite2 if the buffer is empty and the request
        // won't fit.
        return hwrite2(fp, buffer, nbytes, 0);
    }

    if (n > nbytes) n = nbytes;
    memcpy(fp->begin, buffer, n);
    fp->begin += n;
    return (n==nbytes)? (ssize_t) n : hwrite2(fp, buffer, nbytes, n);
}

/// For writing streams, flush buffered output to the underlying stream
/** @return  0 if successful, or `EOF` if an error occurred.

This includes low-level flushing such as via `fdatasync(2)`.
*/
int hflush(hFILE *fp) HTS_RESULT_USED;

/// For hfile_mem: get the internal buffer and it's size from a hfile
/** @return  buffer if successful, or NULL if an error occurred

The buffer returned should not be freed as this will happen when the
hFILE is closed.
*/
char *hfile_mem_get_buffer(hFILE *file, size_t *length);

/// For hfile_mem: get the internal buffer and it's size from a hfile.
/** @return  buffer if successful, or NULL if an error occurred

This is similar to hfile_mem_get_buffer except that ownership of the
buffer is granted to the caller, who now has responsibility for freeing
it.  From this point onwards, the hFILE should not be used for any
purpose other than closing.
*/
char *hfile_mem_steal_buffer(hFILE *file, size_t *length);

#ifdef __cplusplus
}
#endif

#endif
