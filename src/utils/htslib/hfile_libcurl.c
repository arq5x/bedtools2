/*  hfile_libcurl.c -- libcurl backend for low-level file streams.

    Copyright (C) 2015-2017 Genome Research Ltd.

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

#include <config.h>

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#ifndef _WIN32
# include <sys/select.h>
#endif

#include "hfile_internal.h"
#ifdef ENABLE_PLUGINS
#include "version.h"
#endif
#include "htslib/hts.h"  // for hts_version() and hts_verbose
#include "htslib/kstring.h"

#include <curl/curl.h>

// Curl-compatible header linked list
typedef struct {
    struct curl_slist *list;
    unsigned int num;
    unsigned int size;
} hdrlist;

typedef struct {
    hdrlist fixed;                   // List of headers supplied at hopen()
    hdrlist extra;                   // List of headers from callback
    hts_httphdr_callback callback;   // Callback to get more headers
    void *callback_data;             // Data to pass to callback
} http_headers;

typedef struct {
    hFILE base;
    CURL *easy;
    CURLM *multi;
    off_t file_size;
    struct {
        union { char *rd; const char *wr; } ptr;
        size_t len;
    } buffer;
    CURLcode final_result;  // easy result code for finished transfers
    // Flags for communicating with libcurl callbacks:
    unsigned paused : 1;    // callback tells us that it has paused transfer
    unsigned closing : 1;   // informs callback that hclose() has been invoked
    unsigned finished : 1;  // wait_perform() tells us transfer is complete
    unsigned perform_again : 1;
    unsigned is_read : 1;   // Opened in read mode
    unsigned can_seek : 1;  // Can (attempt to) seek on this handle
    int nrunning;
    http_headers headers;
} hFILE_libcurl;

static int http_status_errno(int status)
{
    if (status >= 500)
        switch (status) {
        case 501: return ENOSYS;
        case 503: return EBUSY;
        case 504: return ETIMEDOUT;
        default:  return EIO;
        }
    else if (status >= 400)
        switch (status) {
        case 401: return EPERM;
        case 403: return EACCES;
        case 404: return ENOENT;
        case 405: return EROFS;
        case 407: return EPERM;
        case 408: return ETIMEDOUT;
        case 410: return ENOENT;
        default:  return EINVAL;
        }
    else return 0;
}

static int easy_errno(CURL *easy, CURLcode err)
{
    long lval;

    switch (err) {
    case CURLE_OK:
        return 0;

    case CURLE_UNSUPPORTED_PROTOCOL:
    case CURLE_URL_MALFORMAT:
        return EINVAL;

#if LIBCURL_VERSION_NUM >= 0x071505
    case CURLE_NOT_BUILT_IN:
        return ENOSYS;
#endif

    case CURLE_COULDNT_RESOLVE_PROXY:
    case CURLE_COULDNT_RESOLVE_HOST:
    case CURLE_FTP_CANT_GET_HOST:
        return EDESTADDRREQ; // Lookup failure

    case CURLE_COULDNT_CONNECT:
    case CURLE_SEND_ERROR:
    case CURLE_RECV_ERROR:
        if (curl_easy_getinfo(easy, CURLINFO_OS_ERRNO, &lval) == CURLE_OK)
            return lval;
        else
            return ECONNABORTED;

    case CURLE_REMOTE_ACCESS_DENIED:
    case CURLE_LOGIN_DENIED:
    case CURLE_TFTP_PERM:
        return EACCES;

    case CURLE_PARTIAL_FILE:
        return EPIPE;

    case CURLE_HTTP_RETURNED_ERROR:
        if (curl_easy_getinfo(easy, CURLINFO_RESPONSE_CODE, &lval) == CURLE_OK)
            return http_status_errno(lval);
        else
            return EIO;

    case CURLE_OUT_OF_MEMORY:
        return ENOMEM;

    case CURLE_OPERATION_TIMEDOUT:
        return ETIMEDOUT;

    case CURLE_RANGE_ERROR:
        return ESPIPE;

    case CURLE_SSL_CONNECT_ERROR:
        // TODO return SSL error buffer messages
        return ECONNABORTED;

    case CURLE_FILE_COULDNT_READ_FILE:
    case CURLE_TFTP_NOTFOUND:
        return ENOENT;

    case CURLE_TOO_MANY_REDIRECTS:
        return ELOOP;

    case CURLE_FILESIZE_EXCEEDED:
        return EFBIG;

    case CURLE_REMOTE_DISK_FULL:
        return ENOSPC;

    case CURLE_REMOTE_FILE_EXISTS:
        return EEXIST;

    default:
        return EIO;
    }
}

static int multi_errno(CURLMcode errm)
{
    switch (errm) {
    case CURLM_CALL_MULTI_PERFORM:
    case CURLM_OK:
        return 0;

    case CURLM_BAD_HANDLE:
    case CURLM_BAD_EASY_HANDLE:
    case CURLM_BAD_SOCKET:
        return EBADF;

    case CURLM_OUT_OF_MEMORY:
        return ENOMEM;

    default:
        return EIO;
    }
}


static struct {
    kstring_t useragent;
    CURLSH *share;
    pthread_mutex_t lock;
} curl = { { 0, 0, NULL }, NULL, PTHREAD_MUTEX_INITIALIZER };

static void share_lock(CURL *handle, curl_lock_data data,
                       curl_lock_access access, void *userptr) {
    pthread_mutex_lock(&curl.lock);
}

static void share_unlock(CURL *handle, curl_lock_data data, void *userptr) {
    pthread_mutex_unlock(&curl.lock);
}


static void libcurl_exit()
{
    if (curl_share_cleanup(curl.share) == CURLSHE_OK)
        curl.share = NULL;

    free(curl.useragent.s);
    curl.useragent.l = curl.useragent.m = 0; curl.useragent.s = NULL;

    curl_global_cleanup();
}

static int append_header(hdrlist *hdrs, const char *data, int dup) {
    if (hdrs->num == hdrs->size) {
        unsigned int new_sz = hdrs->size ? hdrs->size * 2 : 4, i;
        struct curl_slist *new_list = realloc(hdrs->list,
                                              new_sz * sizeof(*new_list));
        if (!new_list) return -1;
        hdrs->size = new_sz;
        hdrs->list = new_list;
        for (i = 1; i < hdrs->num; i++) hdrs->list[i-1].next = &hdrs->list[i];
    }
    // Annoyingly, libcurl doesn't declare the char * as const...
    hdrs->list[hdrs->num].data = dup ? strdup(data) : (char *) data;
    if (!hdrs->list[hdrs->num].data) return -1;
    if (hdrs->num > 0) hdrs->list[hdrs->num - 1].next = &hdrs->list[hdrs->num];
    hdrs->list[hdrs->num].next = NULL;
    hdrs->num++;
    return 0;
}

static void free_headers(hdrlist *hdrs, int completely) {
    unsigned int i;
    for (i = 0; i < hdrs->num; i++) {
        free(hdrs->list[i].data);
        hdrs->list[i].data = NULL;
        hdrs->list[i].next = NULL;
    }
    hdrs->num = 0;
    if (completely) {
        free(hdrs->list);
        hdrs->size = 0;
        hdrs->list = NULL;
    }
}

static struct curl_slist * get_header_list(hFILE_libcurl *fp) {
    if (fp->headers.fixed.num > 0)
        return &fp->headers.fixed.list[0];
    if (fp->headers.extra.num > 0)
        return &fp->headers.extra.list[0];
    return 0;
}

static int add_callback_headers(hFILE_libcurl *fp) {
    char **hdrs = NULL, **hdr;

    if (!fp->headers.callback)
        return 0;

    // Get the headers from the callback
    if (fp->headers.callback(fp->headers.callback_data, &hdrs) != 0) {
        return -1;
    }

    if (!hdrs) // No change
        return 0;

    // Remove any old callback headers
    if (fp->headers.fixed.num > 0) {
        // Unlink lists
        fp->headers.fixed.list[fp->headers.fixed.num - 1].next = NULL;
    }
    free_headers(&fp->headers.extra, 0);

    // Convert to libcurl-suitable form
    for (hdr = hdrs; *hdr; hdr++) {
        if (append_header(&fp->headers.extra, *hdr, 0) < 0) {
            goto cleanup;
        }
    }
    for (hdr = hdrs; *hdr; hdr++) *hdr = NULL;

    if (fp->headers.fixed.num > 0 && fp->headers.extra.num > 0) {
        // Relink lists
        fp->headers.fixed.list[fp->headers.fixed.num - 1].next
            = &fp->headers.extra.list[0];
    }
    return 0;

 cleanup:
    while (hdr && *hdr) {
        free(*hdr);
        *hdr = NULL;
    }
    return -1;
}

static void process_messages(hFILE_libcurl *fp)
{
    CURLMsg *msg;
    int remaining;

    while ((msg = curl_multi_info_read(fp->multi, &remaining)) != NULL) {
        switch (msg->msg) {
        case CURLMSG_DONE:
            fp->finished = 1;
            fp->final_result = msg->data.result;
            break;

        default:
            break;
        }
    }
}

static int wait_perform(hFILE_libcurl *fp)
{
    fd_set rd, wr, ex;
    int maxfd, nrunning;
    long timeout;
    CURLMcode errm;

    FD_ZERO(&rd);
    FD_ZERO(&wr);
    FD_ZERO(&ex);
    if (curl_multi_fdset(fp->multi, &rd, &wr, &ex, &maxfd) != CURLM_OK)
        maxfd = -1, timeout = 1000;
    else if (maxfd < 0)
        timeout = 100;  // as recommended by curl_multi_fdset(3)
    else {
        if (curl_multi_timeout(fp->multi, &timeout) != CURLM_OK)
            timeout = 1000;
        else if (timeout < 0)
            timeout = 10000;  // as recommended by curl_multi_timeout(3)
    }

    if (timeout > 0 && ! fp->perform_again) {
        struct timeval tval;
        tval.tv_sec  = (timeout / 1000);
        tval.tv_usec = (timeout % 1000) * 1000;

        if (select(maxfd + 1, &rd, &wr, &ex, &tval) < 0) return -1;
    }

    errm = curl_multi_perform(fp->multi, &nrunning);
    fp->perform_again = 0;
    if (errm == CURLM_CALL_MULTI_PERFORM) fp->perform_again = 1;
    else if (errm != CURLM_OK) { errno = multi_errno(errm); return -1; }

    if (nrunning < fp->nrunning) process_messages(fp);
    return 0;
}


static size_t recv_callback(char *ptr, size_t size, size_t nmemb, void *fpv)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    size_t n = size * nmemb;

    if (n > fp->buffer.len) { fp->paused = 1; return CURL_WRITEFUNC_PAUSE; }
    else if (n == 0) return 0;

    memcpy(fp->buffer.ptr.rd, ptr, n);
    fp->buffer.ptr.rd += n;
    fp->buffer.len -= n;
    return n;
}

static ssize_t libcurl_read(hFILE *fpv, void *bufferv, size_t nbytes)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    char *buffer = (char *) bufferv;
    CURLcode err;

    fp->buffer.ptr.rd = buffer;
    fp->buffer.len = nbytes;
    fp->paused = 0;
    err = curl_easy_pause(fp->easy, CURLPAUSE_CONT);
    if (err != CURLE_OK) { errno = easy_errno(fp->easy, err); return -1; }

    while (! fp->paused && ! fp->finished)
        if (wait_perform(fp) < 0) return -1;

    nbytes = fp->buffer.ptr.rd - buffer;
    fp->buffer.ptr.rd = NULL;
    fp->buffer.len = 0;

    if (fp->finished && fp->final_result != CURLE_OK) {
        errno = easy_errno(fp->easy, fp->final_result);
        return -1;
    }

    return nbytes;
}

static size_t send_callback(char *ptr, size_t size, size_t nmemb, void *fpv)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    size_t n = size * nmemb;

    if (fp->buffer.len == 0) {
        // Send buffer is empty; normally pause, or signal EOF if we're closing
        if (fp->closing) return 0;
        else { fp->paused = 1; return CURL_READFUNC_PAUSE; }
    }

    if (n > fp->buffer.len) n = fp->buffer.len;
    memcpy(ptr, fp->buffer.ptr.wr, n);
    fp->buffer.ptr.wr += n;
    fp->buffer.len -= n;
    return n;
}

static ssize_t libcurl_write(hFILE *fpv, const void *bufferv, size_t nbytes)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    const char *buffer = (const char *) bufferv;
    CURLcode err;

    fp->buffer.ptr.wr = buffer;
    fp->buffer.len = nbytes;
    fp->paused = 0;
    err = curl_easy_pause(fp->easy, CURLPAUSE_CONT);
    if (err != CURLE_OK) { errno = easy_errno(fp->easy, err); return -1; }

    while (! fp->paused && ! fp->finished)
        if (wait_perform(fp) < 0) return -1;

    nbytes = fp->buffer.ptr.wr - buffer;
    fp->buffer.ptr.wr = NULL;
    fp->buffer.len = 0;

    if (fp->finished && fp->final_result != CURLE_OK) {
        errno = easy_errno(fp->easy, fp->final_result);
        return -1;
    }

    return nbytes;
}

static off_t libcurl_seek(hFILE *fpv, off_t offset, int whence)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    hFILE_libcurl temp_fp;
    CURLcode err;
    CURLMcode errm;
    off_t origin, pos;

    if (!fp->is_read || !fp->can_seek) {
        // Cowardly refuse to seek when writing or a previous seek failed.
        errno = ESPIPE;
        return -1;
    }

    switch (whence) {
    case SEEK_SET:
        origin = 0;
        break;
    case SEEK_CUR:
        errno = ENOSYS;
        return -1;
    case SEEK_END:
        if (fp->file_size < 0) { errno = ESPIPE; return -1; }
        origin = fp->file_size;
        break;
    default:
        errno = EINVAL;
        return -1;
    }

    // Check 0 <= origin+offset < fp->file_size carefully, avoiding overflow
    if ((offset < 0)? origin + offset < 0
                : (fp->file_size >= 0 && offset > fp->file_size - origin)) {
        errno = EINVAL;
        return -1;
    }

    pos = origin + offset;

    // TODO If we seem to be doing random access, use CURLOPT_RANGE to do
    // limited reads (e.g. about a BAM block!) so seeking can reuse the
    // existing connection more often.

    // Get new headers from the callback (if defined).  This changes the
    // headers in fp before it gets duplicated, but they should be have been
    // sent by now.

    if (fp->headers.callback) {
        struct curl_slist *list;
        if (add_callback_headers(fp) != 0)
            return -1;
        list = get_header_list(fp);
        if (list) {
            err = curl_easy_setopt(fp->easy, CURLOPT_HTTPHEADER, list);
            if (err != CURLE_OK) {
                errno = easy_errno(fp->easy,err);
                return -1;
            }
        }
    }

    /*
      Duplicate the easy handle, and use CURLOPT_RESUME_FROM_LARGE to open
      a new request to the server, reading from the location that we want
      to seek to.  If the new request works and returns the correct data,
      the original easy handle in *fp is closed and replaced with the new
      one.  If not, we close the new handle, leave *fp unchanged, set
      errno to ESPIPE and return -1 so that the caller knows we can't seek.
      This allows the caller to decide if it wants to continue reading from
      fp, in the same way as it would if reading from a pipe.
     */

    memcpy(&temp_fp, fp, sizeof(temp_fp));
    temp_fp.buffer.len = 0;
    temp_fp.buffer.ptr.rd = NULL;
    temp_fp.easy = curl_easy_duphandle(fp->easy);
    if (!temp_fp.easy)
        goto early_error;

    err = curl_easy_setopt(temp_fp.easy, CURLOPT_RESUME_FROM_LARGE,(curl_off_t)pos);
    err |= curl_easy_setopt(temp_fp.easy, CURLOPT_PRIVATE, &temp_fp);
    err |= curl_easy_setopt(temp_fp.easy, CURLOPT_WRITEDATA, &temp_fp);
    if (err != CURLE_OK)
        goto error;

    temp_fp.buffer.len = 0;  // Ensures we only read the response headers
    temp_fp.paused = temp_fp.finished = 0;

    // fp->multi and temp_fp.multi are the same.
    errm = curl_multi_add_handle(fp->multi, temp_fp.easy);
    if (errm != CURLM_OK) { errno = multi_errno(errm); return -1; }
    temp_fp.nrunning = ++fp->nrunning;

    err = curl_easy_pause(temp_fp.easy, CURLPAUSE_CONT);
    if (err != CURLE_OK)
        goto error_remove;

    while (! temp_fp.paused && ! temp_fp.finished)
        if (wait_perform(&temp_fp) < 0) goto error_remove;

    if (temp_fp.finished && temp_fp.final_result != CURLE_OK)
        goto error_remove;

    // We've got a good response, close the original connection and
    // replace it with the new one.

    errm = curl_multi_remove_handle(fp->multi, fp->easy);
    if (errm != CURLM_OK) {
        curl_easy_reset(temp_fp.easy);
        curl_multi_remove_handle(fp->multi, temp_fp.easy);
        errno = multi_errno(errm);
        return -1;
    }
    fp->nrunning--;

    curl_easy_cleanup(fp->easy);
    fp->easy = temp_fp.easy;
    err = curl_easy_setopt(fp->easy, CURLOPT_WRITEDATA, fp);
    err |= curl_easy_setopt(fp->easy, CURLOPT_PRIVATE, fp);
    if (err != CURLE_OK) {
        int save_errno = easy_errno(fp->easy, err);
        curl_easy_reset(fp->easy);
        errno = save_errno;
        return -1;
    }
    fp->buffer.len = 0;
    fp->paused = temp_fp.paused;
    fp->finished = temp_fp.finished;
    fp->perform_again = temp_fp.perform_again;
    fp->final_result = temp_fp.final_result;

    return pos;

 error_remove:
    curl_easy_reset(temp_fp.easy); // Ensure no pointers to on-stack temp_fp
    errm = curl_multi_remove_handle(fp->multi, temp_fp.easy);
    if (errm != CURLM_OK) {
        errno = multi_errno(errm);
        return -1;
    }
    fp->nrunning--;
 error:
    curl_easy_cleanup(temp_fp.easy);
 early_error:
    fp->can_seek = 0;  // Don't try to seek again
    /* This value for errno may not be entirely true, but the caller may be
       able to carry on with the existing handle. */
    errno = ESPIPE;
    return -1;
}

static int libcurl_close(hFILE *fpv)
{
    hFILE_libcurl *fp = (hFILE_libcurl *) fpv;
    CURLcode err;
    CURLMcode errm;
    int save_errno = 0;

    // Before closing the file, unpause it and perform on it so that uploads
    // have the opportunity to signal EOF to the server -- see send_callback().

    fp->buffer.len = 0;
    fp->closing = 1;
    fp->paused = 0;
    err = curl_easy_pause(fp->easy, CURLPAUSE_CONT);
    if (err != CURLE_OK) save_errno = easy_errno(fp->easy, err);

    while (save_errno == 0 && ! fp->paused && ! fp->finished)
        if (wait_perform(fp) < 0) save_errno = errno;

    if (fp->finished && fp->final_result != CURLE_OK)
        save_errno = easy_errno(fp->easy, fp->final_result);

    errm = curl_multi_remove_handle(fp->multi, fp->easy);
    if (errm != CURLM_OK && save_errno == 0) save_errno = multi_errno(errm);
    fp->nrunning--;

    curl_easy_cleanup(fp->easy);
    curl_multi_cleanup(fp->multi);

    if (fp->headers.callback) // Tell callback to free any data it needs to
        fp->headers.callback(fp->headers.callback_data, NULL);
    free_headers(&fp->headers.fixed, 1);
    free_headers(&fp->headers.extra, 1);

    if (save_errno) { errno = save_errno; return -1; }
    else return 0;
}

static const struct hFILE_backend libcurl_backend =
{
    libcurl_read, libcurl_write, libcurl_seek, NULL, libcurl_close
};

static hFILE *
libcurl_open(const char *url, const char *modes, http_headers *headers)
{
    hFILE_libcurl *fp;
    struct curl_slist *list;
    char mode;
    const char *s;
    CURLcode err;
    CURLMcode errm;
    int save;

    if ((s = strpbrk(modes, "rwa+")) != NULL) {
        mode = *s;
        if (strpbrk(&s[1], "rwa+")) mode = 'e';
    }
    else mode = '\0';

    if (mode != 'r' && mode != 'w') { errno = EINVAL; goto early_error; }

    fp = (hFILE_libcurl *) hfile_init(sizeof (hFILE_libcurl), modes, 0);
    if (fp == NULL) goto early_error;

    if (headers) {
        fp->headers = *headers;
    } else {
        memset(&fp->headers, 0, sizeof(fp->headers));
    }
    fp->file_size = -1;
    fp->buffer.ptr.rd = NULL;
    fp->buffer.len = 0;
    fp->final_result = (CURLcode) -1;
    fp->paused = fp->closing = fp->finished = fp->perform_again = 0;
    fp->can_seek = 1;
    fp->nrunning = 0;
    fp->easy = NULL;

    fp->multi = curl_multi_init();
    if (fp->multi == NULL) { errno = ENOMEM; goto error; }

    fp->easy = curl_easy_init();
    if (fp->easy == NULL) { errno = ENOMEM; goto error; }

    // Make a route to the hFILE_libcurl* given just a CURL* easy handle
    err = curl_easy_setopt(fp->easy, CURLOPT_PRIVATE, fp);

    if (mode == 'r') {
        err |= curl_easy_setopt(fp->easy, CURLOPT_WRITEFUNCTION, recv_callback);
        err |= curl_easy_setopt(fp->easy, CURLOPT_WRITEDATA, fp);
        fp->is_read = 1;
    }
    else {
        err |= curl_easy_setopt(fp->easy, CURLOPT_READFUNCTION, send_callback);
        err |= curl_easy_setopt(fp->easy, CURLOPT_READDATA, fp);
        err |= curl_easy_setopt(fp->easy, CURLOPT_UPLOAD, 1L);
        if (append_header(&fp->headers.fixed,
                          "Transfer-Encoding: chunked", 1) < 0)
            goto error;
        fp->is_read = 0;
    }

    err |= curl_easy_setopt(fp->easy, CURLOPT_SHARE, curl.share);
    err |= curl_easy_setopt(fp->easy, CURLOPT_URL, url);
    err |= curl_easy_setopt(fp->easy, CURLOPT_USERAGENT, curl.useragent.s);
    if (fp->headers.callback) {
        if (add_callback_headers(fp) != 0) goto error;
    }
    if ((list = get_header_list(fp)) != NULL)
        err |= curl_easy_setopt(fp->easy, CURLOPT_HTTPHEADER, list);
    err |= curl_easy_setopt(fp->easy, CURLOPT_FOLLOWLOCATION, 1L);
    if (hts_verbose <= 8)
        err |= curl_easy_setopt(fp->easy, CURLOPT_FAILONERROR, 1L);
    if (hts_verbose >= 8)
        err |= curl_easy_setopt(fp->easy, CURLOPT_VERBOSE, 1L);

    if (err != 0) { errno = ENOSYS; goto error; }

    errm = curl_multi_add_handle(fp->multi, fp->easy);
    if (errm != CURLM_OK) { errno = multi_errno(errm); goto error; }
    fp->nrunning++;

    while (! fp->paused && ! fp->finished)
        if (wait_perform(fp) < 0) goto error_remove;

    if (fp->finished && fp->final_result != CURLE_OK) {
        errno = easy_errno(fp->easy, fp->final_result);
        goto error_remove;
    }

    if (mode == 'r') {
        double dval;
        if (curl_easy_getinfo(fp->easy, CURLINFO_CONTENT_LENGTH_DOWNLOAD,
                              &dval) == CURLE_OK && dval >= 0.0)
            fp->file_size = (off_t) (dval + 0.1);
    }

    fp->base.backend = &libcurl_backend;
    return &fp->base;

error_remove:
    save = errno;
    (void) curl_multi_remove_handle(fp->multi, fp->easy);
    fp->nrunning--;
    errno = save;

error:
    save = errno;
    if (fp->easy) curl_easy_cleanup(fp->easy);
    if (fp->multi) curl_multi_cleanup(fp->multi);
    free_headers(&fp->headers.fixed, 1);
    free_headers(&fp->headers.extra, 1);
    hfile_destroy((hFILE *) fp);
    errno = save;
    return NULL;

early_error:
    save = errno;
    errno = save;
    return NULL;
}

static hFILE *hopen_libcurl(const char *url, const char *modes)
{
    return libcurl_open(url, modes, NULL);
}

static int parse_va_list(http_headers *headers, va_list args)
{
    const char *argtype;

    while ((argtype = va_arg(args, const char *)) != NULL)
        if (strcmp(argtype, "httphdr:v") == 0) {
            const char **hdr;
            for (hdr = va_arg(args, const char **); *hdr; hdr++) {
                if (append_header(&headers->fixed, *hdr, 1) < 0)
                    return -1;
            }
        }
        else if (strcmp(argtype, "httphdr:l") == 0) {
            const char *hdr;
            while ((hdr = va_arg(args, const char *)) != NULL) {
                if (append_header(&headers->fixed, hdr, 1) < 0)
                    return -1;
            }
        }
        else if (strcmp(argtype, "httphdr") == 0) {
            const char *hdr = va_arg(args, const char *);
            if (hdr) {
                if (append_header(&headers->fixed, hdr, 1) < 0)
                    return -1;
            }
        }
        else if (strcmp(argtype, "httphdr_callback") == 0) {
            headers->callback = va_arg(args, const hts_httphdr_callback);
        }
        else if (strcmp(argtype, "httphdr_callback_data") == 0) {
            headers->callback_data = va_arg(args, void *);
        }
        else if (strcmp(argtype, "va_list") == 0) {
            va_list *args2 = va_arg(args, va_list *);
            if (args2) {
                if (parse_va_list(headers, *args2) < 0) return -1;
            }
        }
        else { errno = EINVAL; return -1; }

    return 0;
}

/*
  HTTP headers to be added to the request can be passed in as extra
  arguments to hopen().  The headers can be specified as follows:

  * Single header:
    hopen(url, mode, "httphdr", "X-Hdr-1: text", NULL);

  * Multiple headers in the argument list:
    hopen(url, mode, "httphdr:l", "X-Hdr-1: text", "X-Hdr-2: text", NULL, NULL);

  * Multiple headers in a char* array:
    hopen(url, mode, "httphdr:v", hdrs, NULL);
    where `hdrs` is a char **.  The list ends with a NULL pointer.

  * A callback function
    hopen(url, mode, "httphdr_callback", func,
                     "httphdr_callback_data", arg, NULL);
    `func` has type
         int (* hts_httphdr_callback) (void *cb_data, char ***hdrs);
    `arg` is passed to the callback as a void *.

    The function is called at file open, and when attempting to seek (which
    opens a new HTTP request).  This allows, for example, access tokens
    that may have gone stale to be regenerated.  The function is also
    called (with `hdrs` == NULL) on file close so that the callback can
    free any memory that it needs to.

    The callback should return 0 on success, non-zero on failure.  It should
    return in *hdrs a list of strings containing the new headers (terminated
    with a NULL pointer).  These will replace any headers previously supplied
    by the callback.  If no changes are necessary, it can return NULL
    in *hdrs, in which case the previous headers will be left unchanged.

    Ownership of the strings in the header list passes to hfile_libcurl,
    so the callback should not attempt to use or free them itself.  The memory
    containing the array belongs to the callback and will not be freed by
    hfile_libcurl.

    Headers supplied by the callback are appended after any specified
    using the "httphdr", "httphdr:l" or "httphdr:v" methods.  No attempt
    is made to replace these headers (even if a key is repeated) so anything
    that is expected to vary needs to come from the callback.
 */

static hFILE *vhopen_libcurl(const char *url, const char *modes, va_list args)
{
    hFILE *fp = NULL;
    http_headers headers = { { NULL, 0, 0 }, { NULL, 0, 0 }, NULL, NULL };
    if (parse_va_list(&headers, args) == 0) {
        fp = libcurl_open(url, modes, &headers);
    }

    if (!fp) {
        free_headers(&headers.fixed, 1);
    }
    return fp;
}

int PLUGIN_GLOBAL(hfile_plugin_init,_libcurl)(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { hopen_libcurl, hfile_always_remote, "libcurl",
          2000 + 50,
          vhopen_libcurl };

#ifdef ENABLE_PLUGINS
    // Embed version string for examination via strings(1) or what(1)
    static const char id[] = "@(#)hfile_libcurl plugin (htslib)\t" HTS_VERSION;
    const char *version = strchr(id, '\t')+1;
#else
    const char *version = hts_version();
#endif
    const curl_version_info_data *info;
    const char * const *protocol;
    CURLcode err;
    CURLSHcode errsh;

    err = curl_global_init(CURL_GLOBAL_ALL);
    if (err != CURLE_OK) { errno = easy_errno(NULL, err); return -1; }

    curl.share = curl_share_init();
    if (curl.share == NULL) { curl_global_cleanup(); errno = EIO; return -1; }
    errsh = curl_share_setopt(curl.share, CURLSHOPT_LOCKFUNC, share_lock);
    errsh |= curl_share_setopt(curl.share, CURLSHOPT_UNLOCKFUNC, share_unlock);
    errsh |= curl_share_setopt(curl.share, CURLSHOPT_SHARE, CURL_LOCK_DATA_DNS);
    if (errsh != 0) {
        curl_share_cleanup(curl.share);
        curl_global_cleanup();
        errno = EIO;
        return -1;
    }

    info = curl_version_info(CURLVERSION_NOW);
    ksprintf(&curl.useragent, "htslib/%s libcurl/%s", version, info->version);

    self->name = "libcurl";
    self->destroy = libcurl_exit;

    for (protocol = info->protocols; *protocol; protocol++)
        hfile_add_scheme_handler(*protocol, &handler);
    return 0;
}
