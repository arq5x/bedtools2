/*  hfile_internal.h -- internal parts of low-level input/output streams.

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

#ifndef HFILE_INTERNAL_H
#define HFILE_INTERNAL_H

#include <stdarg.h>

#include "htslib/hfile.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
  @abstract  Resizes the buffer within an hFILE.

  @notes  Changes the buffer size for an hFILE.  Ideally this is done
  immediately after opening.  If performed later, this function may
  fail if we are reducing the buffer size and the current offset into
  the buffer is beyond the new capacity.

  @param fp        The file stream
  @param bufsiz    The size of the new buffer

  @return Returns 0 on success, -1 on failure.
 */
int hfile_set_blksize(hFILE *fp, size_t bufsiz);

struct BGZF;
/*!
  @abstract Return the hFILE connected to a BGZF
 */
struct hFILE *bgzf_hfile(struct BGZF *fp);

struct hFILE_backend {
    /* As per read(2), returning the number of bytes read (possibly 0) or
       negative (and setting errno) on errors.  Front-end code will call this
       repeatedly if necessary to attempt to get the desired byte count.  */
    ssize_t (*read)(hFILE *fp, void *buffer, size_t nbytes) HTS_RESULT_USED;

    /* As per write(2), returning the number of bytes written or negative (and
       setting errno) on errors.  Front-end code will call this repeatedly if
       necessary until the desired block is written or an error occurs.  */
    ssize_t (*write)(hFILE *fp, const void *buffer, size_t nbytes)
        HTS_RESULT_USED;

    /* As per lseek(2), returning the resulting offset within the stream or
       negative (and setting errno) on errors.  */
    off_t (*seek)(hFILE *fp, off_t offset, int whence) HTS_RESULT_USED;

    /* Performs low-level flushing, if any, e.g., fsync(2); for writing streams
       only.  Returns 0 for success or negative (and sets errno) on errors. */
    int (*flush)(hFILE *fp) HTS_RESULT_USED;

    /* Closes the underlying stream (for output streams, the buffer will
       already have been flushed), returning 0 for success or negative (and
       setting errno) on errors, as per close(2).  */
    int (*close)(hFILE *fp) HTS_RESULT_USED;
};

/* May be called by hopen_*() functions to decode a fopen()-style mode into
   open(2)-style flags.  */
int hfile_oflags(const char *mode);

/* Must be called by hopen_*() functions to allocate the hFILE struct and set
   up its base.  Capacity is a suggested buffer size (e.g., via fstat(2))
   or 0 for a default-sized buffer.  */
hFILE *hfile_init(size_t struct_size, const char *mode, size_t capacity);

/* Alternative to hfile_init() for in-memory backends for which the base
   buffer is the only storage.  Buffer is already allocated via malloc(2)
   of size buf_size and with buf_filled bytes already filled.  Ownership
   of the buffer is transferred to the resulting hFILE.  */
hFILE *hfile_init_fixed(size_t struct_size, const char *mode,
                        char *buffer, size_t buf_filled, size_t buf_size);

/* May be called by hopen_*() functions to undo the effects of hfile_init()
   in the event opening the stream subsequently fails.  (This is safe to use
   even if fp is NULL.  This takes care to preserve errno.)  */
void hfile_destroy(hFILE *fp);


struct hFILE_scheme_handler {
    /* Opens a stream when dispatched by hopen(); should call hfile_init()
       to malloc a struct "derived" from hFILE and initialise it appropriately,
       including setting base.backend to its own backend vector.  */
    hFILE *(*open)(const char *filename, const char *mode) HTS_RESULT_USED;

    /* Returns whether the URL denotes remote storage when dispatched by
       hisremote().  For simple cases, use one of hfile_always_*() below.  */
    int (*isremote)(const char *filename) HTS_RESULT_USED;

    /* The name of the plugin or other code providing this handler.  */
    const char *provider;

    /* If multiple handlers are registered for the same scheme, the one with
       the highest priority is used; range is 0 (lowest) to 100 (highest).
       This field is used modulo 1000 as a priority; thousands indicate
       later revisions to this structure, as noted below.  */
    int priority;

    /* Fields below are present when priority >= 2000.  */

    /* Same as the open() method, used when extra arguments have been given
       to hopen().  */
    hFILE *(*vopen)(const char *filename, const char *mode, va_list args)
        HTS_RESULT_USED;
};

/* May be used as an isremote() function in simple cases.  */
extern int hfile_always_local (const char *fname);
extern int hfile_always_remote(const char *fname);

/* Should be called by plugins for each URL scheme they wish to handle.  */
void hfile_add_scheme_handler(const char *scheme,
                              const struct hFILE_scheme_handler *handler);

struct hFILE_plugin {
    /* On entry, HTSlib's plugin API version (currently 1).  */
    int api_version;

    /* On entry, the plugin's handle as returned by dlopen() etc.  */
    void *obj;

    /* The plugin should fill this in with its (human-readable) name.  */
    const char *name;

    /* The plugin may wish to fill in a function to be called on closing.  */
    void (*destroy)(void);
};

#ifdef ENABLE_PLUGINS
#define PLUGIN_GLOBAL(identifier,suffix) identifier

/* Plugins must define an entry point with this signature.  */
extern int hfile_plugin_init(struct hFILE_plugin *self);

#else
#define PLUGIN_GLOBAL(identifier,suffix) identifier##suffix

/* Only plugins distributed within the HTSlib source that might be built
   even with --disable-plugins need to use PLUGIN_GLOBAL and be listed here;
   others can simply define hfile_plugin_init().  */

extern int hfile_plugin_init_gcs(struct hFILE_plugin *self);
extern int hfile_plugin_init_libcurl(struct hFILE_plugin *self);
extern int hfile_plugin_init_s3(struct hFILE_plugin *self);
#endif

/* This one is never built as a separate plugin.  */
extern int hfile_plugin_init_net(struct hFILE_plugin *self);

// Callback to allow headers to be set in http connections.  Currently used
// to allow s3 to renew tokens when seeking.  Kept internal for now,
// although we may consider exposing it in the API later.
typedef int (* hts_httphdr_callback) (void *cb_data, char ***hdrs);

#ifdef __cplusplus
}
#endif

#endif
