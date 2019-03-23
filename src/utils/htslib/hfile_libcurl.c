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
#include <assert.h>

#include "hfile_internal.h"
#ifdef ENABLE_PLUGINS
#include "version.h"
#endif
#include "htslib/hts.h"  // for hts_version() and hts_verbose
#include "htslib/kstring.h"
#include "htslib/khash.h"

#include <curl/curl.h>

// Number of seconds to take off auth_token expiry, to allow for clock skew
// and slow servers
#define AUTH_REFRESH_EARLY_SECS 60

// Minimum number of bytes to skip when seeking forward.  Seeks less than
// this will just read the data and throw it away.  The optimal value
// depends on how long it takes to make a new connection compared
// to how fast the data arrives.
#define MIN_SEEK_FORWARD 1000000

typedef struct {
    char *path;
    char *token;
    time_t expiry;
    int failed;
    pthread_mutex_t lock;
} auth_token;

// For the authorization header cache
KHASH_MAP_INIT_STR(auth_map, auth_token *)

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
    auth_token *auth;                // Authentication token
    int auth_hdr_num;                // Location of auth_token in hdrlist extra
                                     // If -1, Authorization header is in fixed
                                     //    -2, it came from the callback
                                     //    -3, "auth_token_enabled", "false"
                                     //        passed to hopen()
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
    unsigned is_recursive:1; // Opened by hfile_libcurl itself
    unsigned tried_seek : 1; // At least one seek has been attempted
    int nrunning;
    http_headers headers;
    off_t delayed_seek;      // Location to seek to before reading
    off_t last_offset;       // Location we're seeking from
} hFILE_libcurl;

static off_t libcurl_seek(hFILE *fpv, off_t offset, int whence);
static int restart_from_position(hFILE_libcurl *fp, off_t pos);

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
    char *auth_path;
    khash_t(auth_map) *auth_map;
    int allow_unencrypted_auth_header;
    pthread_mutex_t auth_lock;
    pthread_mutex_t share_lock;
} curl = { { 0, 0, NULL }, NULL, NULL, NULL, 0, PTHREAD_MUTEX_INITIALIZER,
           PTHREAD_MUTEX_INITIALIZER };

static void share_lock(CURL *handle, curl_lock_data data,
                       curl_lock_access access, void *userptr) {
    pthread_mutex_lock(&curl.share_lock);
}

static void share_unlock(CURL *handle, curl_lock_data data, void *userptr) {
    pthread_mutex_unlock(&curl.share_lock);
}

static void free_auth(auth_token *tok) {
    if (!tok) return;
    if (pthread_mutex_destroy(&tok->lock)) abort();
    free(tok->path);
    free(tok->token);
    free(tok);
}

static void libcurl_exit()
{
    if (curl_share_cleanup(curl.share) == CURLSHE_OK)
        curl.share = NULL;

    free(curl.useragent.s);
    curl.useragent.l = curl.useragent.m = 0; curl.useragent.s = NULL;

    free(curl.auth_path);
    curl.auth_path = NULL;

    if (curl.auth_map) {
        khiter_t i;
        for (i = kh_begin(curl.auth_map); i != kh_end(curl.auth_map); ++i) {
            if (kh_exist(curl.auth_map, i)) {
                free_auth(kh_value(curl.auth_map, i));
                kh_key(curl.auth_map, i) = NULL;
                kh_value(curl.auth_map, i) = NULL;
            }
        }
        kh_destroy(auth_map, curl.auth_map);
        curl.auth_map = NULL;
    }

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

static inline int is_authorization(const char *hdr) {
    return (strncasecmp("authorization:", hdr, 14) == 0);
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

    if (fp->headers.auth_hdr_num > 0 || fp->headers.auth_hdr_num == -2)
        fp->headers.auth_hdr_num = 0; // Just removed it...

    // Convert to libcurl-suitable form
    for (hdr = hdrs; *hdr; hdr++) {
        if (append_header(&fp->headers.extra, *hdr, 0) < 0) {
            goto cleanup;
        }
        if (is_authorization(*hdr) && !fp->headers.auth_hdr_num)
            fp->headers.auth_hdr_num = -2;
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

/*
 * Read an OAUTH2-style Bearer access token (see
 * https://tools.ietf.org/html/rfc6750#section-4).
 * Returns 'v' for valid; 'i' for invalid (token missing or wrong sort);
 * '?' for a JSON parse error; 'm' if it runs out of memory.
 */
static int read_auth_json(auth_token *tok, hFILE *auth_fp) {
    hts_json_token *t = hts_json_alloc_token();
    kstring_t str = {0, 0, NULL};
    char *token = NULL, *type = NULL, *expiry = NULL;
    int ret = 'i';

    if (!t) goto error;

    if ((ret = hts_json_fnext(auth_fp, t, &str)) != '{') goto error;
    while (hts_json_fnext(auth_fp, t, &str) != '}') {
        char *key;
        if (hts_json_token_type(t) != 's') {
            ret = '?';
            goto error;
        }
        key = hts_json_token_str(t);
        if (!key) goto error;
        if (strcmp(key, "access_token") == 0) {
            if ((ret = hts_json_fnext(auth_fp, t, &str)) != 's') goto error;
            token = ks_release(&str);
        } else if (strcmp(key, "token_type") == 0) {
            if ((ret = hts_json_fnext(auth_fp, t, &str)) != 's') goto error;
            type = ks_release(&str);
        } else if (strcmp(key, "expires_in") == 0) {
            if ((ret = hts_json_fnext(auth_fp, t, &str)) != 'n') goto error;
            expiry = ks_release(&str);
        } else if (hts_json_fskip_value(auth_fp, '\0') != 'v') {
            ret = '?';
            goto error;
        }
    }

    if (!token || (type && strcmp(type, "Bearer") != 0)) {
        ret = 'i';
        goto error;
    }

    ret = 'm';
    str.l = 0;
    if (kputs("Authorization: Bearer ", &str) < 0) goto error;
    if (kputs(token, &str) < 0) goto error;
    free(tok->token);
    tok->token = ks_release(&str);
    if (expiry) {
        long exp = strtol(expiry, NULL, 10);
        if (exp < 0) exp = 0;
        tok->expiry = time(NULL) + exp;
    } else {
        tok->expiry = 0;
    }
    ret = 'v';

 error:
    free(token);
    free(type);
    free(expiry);
    free(str.s);
    hts_json_free_token(t);
    return ret;
}

static int read_auth_plain(auth_token *tok, hFILE *auth_fp) {
    kstring_t line = {0, 0, NULL};
    kstring_t token = {0, 0, NULL};
    const char *start, *end;

    if (kgetline(&line, (char * (*)(char *, int, void *)) hgets, auth_fp) < 0) goto error;
    if (kputc('\0', &line) < 0) goto error;

    for (start = line.s; *start && isspace_c(*start); start++) {}
    for (end = start; *end && !isspace_c(*end); end++) {}

    if (end > start) {
        if (kputs("Authorization: Bearer ", &token) < 0) goto error;
        if (kputsn(start, end - start, &token) < 0) goto error;
    }

    free(tok->token);
    tok->token = ks_release(&token);
    tok->expiry = 0;
    free(line.s);
    return 0;

 error:
    free(line.s);
    free(token.s);
    return -1;
}

static int renew_auth_token(auth_token *tok, int *changed) {
    hFILE *auth_fp = NULL;
    char buffer[16];
    ssize_t len;

    *changed = 0;
    if (tok->expiry == 0 || time(NULL) + AUTH_REFRESH_EARLY_SECS < tok->expiry)
        return 0; // Still valid

    if (tok->failed)
        return -1;

    *changed = 1;
    auth_fp = hopen(tok->path, "rR");
    if (!auth_fp) {
        // Not worried about missing files; other errors are bad.
        if (errno != ENOENT)
            goto fail;

        tok->expiry = 0; // Prevent retry
        free(tok->token); // Just in case it was set
        return 0;
    }

    len = hpeek(auth_fp, buffer, sizeof(buffer));
    if (len < 0)
        goto fail;

    if (memchr(buffer, '{', len) != NULL) {
        if (read_auth_json(tok, auth_fp) != 'v')
            goto fail;
    } else {
        if (read_auth_plain(tok, auth_fp) < 0)
            goto fail;
    }

    return hclose(auth_fp) < 0 ? -1 : 0;

 fail:
    tok->failed = 1;
    if (auth_fp) hclose_abruptly(auth_fp);
    return -1;
}

static int add_auth_header(hFILE_libcurl *fp) {
    int changed = 0;

    if (fp->headers.auth_hdr_num < 0)
        return 0; // Have an Authorization header from open or header callback

    if (!fp->headers.auth)
        return 0; // Nothing to add

    pthread_mutex_lock(&fp->headers.auth->lock);
    if (renew_auth_token(fp->headers.auth, &changed) < 0)
        goto unlock_fail;

    if (!changed && fp->headers.auth_hdr_num > 0) {
        pthread_mutex_unlock(&fp->headers.auth->lock);
        return 0;
    }

    if (fp->headers.auth_hdr_num > 0) {
        // Had a previous header, so swap in the new one
        char *header = fp->headers.auth->token;
        char *header_copy = header ? strdup(header) : NULL;
        int idx = fp->headers.auth_hdr_num - 1;
        if (header && !header_copy)
            goto unlock_fail;

        if (header_copy) {
            free(fp->headers.extra.list[idx].data);
            fp->headers.extra.list[idx].data = header_copy;
        } else {
            unsigned int j;
            // More complicated case - need to get rid of the old header
            // and tidy up linked lists
            free(fp->headers.extra.list[idx].data);
            for (j = idx + 1; j < fp->headers.extra.num; j++) {
                fp->headers.extra.list[j - 1] = fp->headers.extra.list[j];
                fp->headers.extra.list[j - 1].next = &fp->headers.extra.list[j];
            }
            fp->headers.extra.num--;
            if (fp->headers.extra.num > 0) {
                fp->headers.extra.list[fp->headers.extra.num-1].next = NULL;
            } else if (fp->headers.fixed.num > 0) {
                fp->headers.fixed.list[fp->headers.fixed.num - 1].next = NULL;
            }
            fp->headers.auth_hdr_num = 0;
        }
    } else if (fp->headers.auth->token) {
        // Add new header and remember where it is
        if (append_header(&fp->headers.extra,
                          fp->headers.auth->token, 1) < 0) {
            goto unlock_fail;
        }
        fp->headers.auth_hdr_num = fp->headers.extra.num;
    }

    pthread_mutex_unlock(&fp->headers.auth->lock);
    return 0;

 unlock_fail:
    pthread_mutex_unlock(&fp->headers.auth->lock);
    return -1;
}

static int get_auth_token(hFILE_libcurl *fp, const char *url) {
    const char *host = NULL, *p, *q;
    kstring_t name = {0, 0, NULL};
    size_t host_len = 0;
    khiter_t idx;
    auth_token *tok = NULL;

    // Nothing to do if:
    //   curl.auth_path has not been set
    //   fp was made by hfile_libcurl (e.g. auth_path is a http:// url)
    //   we already have an Authorization header
    if (!curl.auth_path || fp->is_recursive || fp->headers.auth_hdr_num != 0)
        return 0;

    // Insist on having a secure connection unless the user insists harder
    if (!curl.allow_unencrypted_auth_header && strncmp(url, "https://", 8) != 0)
        return 0;

    host = strstr(url, "://");
    if (host) {
        host += 3;
        host_len = strcspn(host, "/");
    }

    p = curl.auth_path;
    while ((q = strstr(p, "%h")) != NULL) {
        if (q - p > INT_MAX || host_len > INT_MAX) goto error;
        if (kputsn_(p, q - p, &name) < 0) goto error;
        if (kputsn_(host, host_len, &name) < 0) goto error;
        p = q + 2;
    }
    if (kputs(p, &name) < 0) goto error;

    pthread_mutex_lock(&curl.auth_lock);
    idx = kh_get(auth_map, curl.auth_map, name.s);
    if (idx < kh_end(curl.auth_map)) {
        tok = kh_value(curl.auth_map, idx);
    } else {
        tok = calloc(1, sizeof(*tok));
        if (tok && pthread_mutex_init(&tok->lock, NULL) != 0) {
            free(tok);
            tok = NULL;
        }
        if (tok) {
            int ret = -1;
            tok->path = ks_release(&name);
            tok->token = NULL;
            tok->expiry = 1; // Force refresh
            idx = kh_put(auth_map, curl.auth_map, tok->path, &ret);
            if (ret < 0) {
                free_auth(tok);
                tok = NULL;
            }
            kh_value(curl.auth_map, idx) = tok;
        }
    }
    pthread_mutex_unlock(&curl.auth_lock);

    fp->headers.auth = tok;
    free(name.s);

    return add_auth_header(fp);

 error:
    free(name.s);
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

    if (!fp->perform_again) {
        FD_ZERO(&rd);
        FD_ZERO(&wr);
        FD_ZERO(&ex);
        if (curl_multi_fdset(fp->multi, &rd, &wr, &ex, &maxfd) != CURLM_OK)
            maxfd = -1, timeout = 1000;
        else {
            if (curl_multi_timeout(fp->multi, &timeout) != CURLM_OK)
                timeout = 1000;
            else if (timeout < 0) {
                timeout = 10000;  // as recommended by curl_multi_timeout(3)
            }
        }
        if (maxfd < 0 && timeout > 100)
            timeout = 100; // as recommended by curl_multi_fdset(3)

        if (timeout > 0) {
            struct timeval tval;
            tval.tv_sec  = (timeout / 1000);
            tval.tv_usec = (timeout % 1000) * 1000;

            if (select(maxfd + 1, &rd, &wr, &ex, &tval) < 0) return -1;
        }
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
    off_t to_skip = -1;
    ssize_t got = 0;
    CURLcode err;

    if (fp->delayed_seek >= 0) {
        assert(fp->base.offset == fp->delayed_seek
               && fp->base.begin == fp->base.buffer
               && fp->base.end == fp->base.buffer);

        if (fp->last_offset >= 0
            && fp->delayed_seek > fp->last_offset
            && fp->delayed_seek - fp->last_offset < MIN_SEEK_FORWARD) {
            // If not seeking far, just read the data and throw it away.  This
            // is likely to be quicker than opening a new stream
            to_skip = fp->delayed_seek - fp->last_offset;
        } else {
            if (restart_from_position(fp, fp->delayed_seek) < 0) {
                return -1;
            }
        }
        fp->delayed_seek = -1;
        fp->last_offset = -1;
    }

    do {
        fp->buffer.ptr.rd = buffer;
        fp->buffer.len = nbytes;
        fp->paused = 0;
        err = curl_easy_pause(fp->easy, CURLPAUSE_CONT);
        if (err != CURLE_OK) { errno = easy_errno(fp->easy, err); return -1; }

        while (! fp->paused && ! fp->finished)
            if (wait_perform(fp) < 0) return -1;

        got = fp->buffer.ptr.rd - buffer;

        if (to_skip >= 0) { // Skipping over a small seek
            if (got < to_skip) { // Need to skip more data
                to_skip -= got;
            } else {
                got -= to_skip;
                if (got > 0) {  // If enough was skipped, return the rest
                    memmove(buffer, buffer + to_skip, got);
                    to_skip = -1;
                }
            }
        }
    } while (to_skip >= 0 && ! fp->finished);
    fp->buffer.ptr.rd = NULL;
    fp->buffer.len = 0;

    if (fp->finished && fp->final_result != CURLE_OK) {
        errno = easy_errno(fp->easy, fp->final_result);
        return -1;
    }

    return got;
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

    if (fp->tried_seek) {
        /* Seeking has worked at least once, so now we can delay doing
           the actual work until the next read.  This avoids lots of pointless
           http or ftp reconnections if the caller does lots of seeks
           without any intervening reads. */
        if (fp->delayed_seek < 0) {
            fp->last_offset = fp->base.offset + (fp->base.end - fp->base.buffer);
        }
        fp->delayed_seek = pos;
        return pos;
    }

    if (restart_from_position(fp, pos) < 0) {
        /* This value for errno may not be entirely true, but the caller may be
           able to carry on with the existing handle. */
        errno = ESPIPE;
        return -1;
    }

    fp->tried_seek = 1;
    return pos;
}

static int restart_from_position(hFILE_libcurl *fp, off_t pos) {
    hFILE_libcurl temp_fp;
    CURLcode err;
    CURLMcode errm;
    int update_headers = 0;
    int save_errno = 0;

    // TODO If we seem to be doing random access, use CURLOPT_RANGE to do
    // limited reads (e.g. about a BAM block!) so seeking can reuse the
    // existing connection more often.

    // Get new headers from the callback (if defined).  This changes the
    // headers in fp before it gets duplicated, but they should be have been
    // sent by now.

    if (fp->headers.callback) {
        if (add_callback_headers(fp) != 0)
            return -1;
        update_headers = 1;
    }
    if (fp->headers.auth_hdr_num > 0 && fp->headers.auth) {
        if (add_auth_header(fp) != 0)
            return -1;
        update_headers = 1;
    }
    if (update_headers) {
        struct curl_slist *list = get_header_list(fp);
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
      one.  If not, we close the new handle and leave *fp unchanged.
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
    if (err != CURLE_OK) {
        save_errno = easy_errno(temp_fp.easy, err);
        goto error;
    }

    temp_fp.buffer.len = 0;  // Ensures we only read the response headers
    temp_fp.paused = temp_fp.finished = 0;

    // fp->multi and temp_fp.multi are the same.
    errm = curl_multi_add_handle(fp->multi, temp_fp.easy);
    if (errm != CURLM_OK) {
        save_errno = multi_errno(errm);
        goto error;
    }
    temp_fp.nrunning = ++fp->nrunning;

    err = curl_easy_pause(temp_fp.easy, CURLPAUSE_CONT);
    if (err != CURLE_OK) {
        save_errno = easy_errno(temp_fp.easy, err);
        goto error_remove;
    }

    while (! temp_fp.paused && ! temp_fp.finished)
        if (wait_perform(&temp_fp) < 0) {
            save_errno = errno;
            goto error_remove;
        }

    if (temp_fp.finished && temp_fp.final_result != CURLE_OK) {
        save_errno = easy_errno(temp_fp.easy, temp_fp.final_result);
        goto error_remove;
    }

    // We've got a good response, close the original connection and
    // replace it with the new one.

    errm = curl_multi_remove_handle(fp->multi, fp->easy);
    if (errm != CURLM_OK) {
        // Clean up as much as possible
        curl_easy_reset(temp_fp.easy);
        if (curl_multi_remove_handle(fp->multi, temp_fp.easy) == CURLM_OK) {
            fp->nrunning--;
            curl_easy_cleanup(temp_fp.easy);
        }
        save_errno = multi_errno(errm);
        goto early_error;
    }
    fp->nrunning--;

    curl_easy_cleanup(fp->easy);
    fp->easy = temp_fp.easy;
    err = curl_easy_setopt(fp->easy, CURLOPT_WRITEDATA, fp);
    err |= curl_easy_setopt(fp->easy, CURLOPT_PRIVATE, fp);
    if (err != CURLE_OK) {
        save_errno = easy_errno(fp->easy, err);
        curl_easy_reset(fp->easy);
        errno = save_errno;
        return -1;
    }
    fp->buffer.len = 0;
    fp->paused = temp_fp.paused;
    fp->finished = temp_fp.finished;
    fp->perform_again = temp_fp.perform_again;
    fp->final_result = temp_fp.final_result;

    return 0;

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
    if (save_errno)
        errno = save_errno;
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
    int save, is_recursive;

    is_recursive = strchr(modes, 'R') != NULL;

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
    fp->tried_seek = 0;
    fp->delayed_seek = fp->last_offset = -1;
    fp->is_recursive = is_recursive;
    fp->nrunning = 0;
    fp->easy = NULL;

    fp->multi = curl_multi_init();
    if (fp->multi == NULL) { errno = ENOMEM; goto error; }

    fp->easy = curl_easy_init();
    if (fp->easy == NULL) { errno = ENOMEM; goto error; }

    // Make a route to the hFILE_libcurl* given just a CURL* easy handle
    err = curl_easy_setopt(fp->easy, CURLOPT_PRIVATE, fp);

    // Avoid many repeated CWD calls with FTP, instead requesting the filename
    // by full path (as done in knet, but not strictly compliant with RFC1738).
    err |= curl_easy_setopt(fp->easy, CURLOPT_FTP_FILEMETHOD, CURLFTPMETHOD_NOCWD);

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
    {
        char* env_curl_ca_bundle = getenv("CURL_CA_BUNDLE");
        if (env_curl_ca_bundle) {
            err |= curl_easy_setopt(fp->easy, CURLOPT_CAINFO, env_curl_ca_bundle);
        }
    }
    err |= curl_easy_setopt(fp->easy, CURLOPT_USERAGENT, curl.useragent.s);
    if (fp->headers.callback) {
        if (add_callback_headers(fp) != 0) goto error;
    }
    if (get_auth_token(fp, url) < 0)
        goto error;
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
    free_headers(&fp->headers.extra, 1);
    hfile_destroy((hFILE *) fp);
    errno = save;
    return NULL;

early_error:
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
                if (is_authorization(*hdr))
                    headers->auth_hdr_num = -1;
            }
        }
        else if (strcmp(argtype, "httphdr:l") == 0) {
            const char *hdr;
            while ((hdr = va_arg(args, const char *)) != NULL) {
                if (append_header(&headers->fixed, hdr, 1) < 0)
                    return -1;
                if (is_authorization(hdr))
                    headers->auth_hdr_num = -1;
            }
        }
        else if (strcmp(argtype, "httphdr") == 0) {
            const char *hdr = va_arg(args, const char *);
            if (hdr) {
                if (append_header(&headers->fixed, hdr, 1) < 0)
                    return -1;
                if (is_authorization(hdr))
                    headers->auth_hdr_num = -1;
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
        else if (strcmp(argtype, "auth_token_enabled") == 0) {
            const char *flag = va_arg(args, const char *);
            if (strcmp(flag, "false") == 0)
                headers->auth_hdr_num = -3;
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
    const char *auth;
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

    if ((auth = getenv("HTS_AUTH_LOCATION")) != NULL) {
        curl.auth_path = strdup(auth);
        curl.auth_map = kh_init(auth_map);
        if (!curl.auth_path || !curl.auth_map) {
            int save_errno = errno;
            free(curl.auth_path);
            kh_destroy(auth_map, curl.auth_map);
            curl_share_cleanup(curl.share);
            curl_global_cleanup();
            errno = save_errno;
            return -1;
        }
    }
    if ((auth = getenv("HTS_ALLOW_UNENCRYPTED_AUTHORIZATION_HEADER")) != NULL
        && strcmp(auth, "I understand the risks") == 0) {
        curl.allow_unencrypted_auth_header = 1;
    }

    info = curl_version_info(CURLVERSION_NOW);
    ksprintf(&curl.useragent, "htslib/%s libcurl/%s", version, info->version);

    self->name = "libcurl";
    self->destroy = libcurl_exit;

    for (protocol = info->protocols; *protocol; protocol++)
        hfile_add_scheme_handler(*protocol, &handler);
    return 0;
}
