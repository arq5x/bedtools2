/*  multipart.c -- GA4GH redirection and multipart backend for file streams.

    Copyright (C) 2016 Genome Research Ltd.

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

#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "htslib/kstring.h"

#include "hts_internal.h"
#include "hfile_internal.h"

#ifndef EPROTO
#define EPROTO ENOEXEC
#endif

typedef struct hfile_part {
    char *url;
    char **headers;
} hfile_part;

typedef struct {
    hFILE base;
    hfile_part *parts;
    size_t nparts, maxparts, current;
    hFILE *currentfp;
} hFILE_multipart;

static void free_part(hfile_part *p)
{
    free(p->url);
    if (p->headers) {
        char **hdr;
        for (hdr = p->headers; *hdr; hdr++) free(*hdr);
        free(p->headers);
    }

    p->url = NULL;
    p->headers = NULL;
}

static void free_all_parts(hFILE_multipart *fp)
{
    size_t i;
    for (i = 0; i < fp->nparts; i++) free_part(&fp->parts[i]);
    free(fp->parts);
}

static ssize_t multipart_read(hFILE *fpv, void *buffer, size_t nbytes)
{
    hFILE_multipart *fp = (hFILE_multipart *) fpv;
    size_t n;

open_next:
    if (fp->currentfp == NULL) {
        if (fp->current < fp->nparts) {
            const hfile_part *p = &fp->parts[fp->current];
            hts_log_debug("Opening part #%zu of %zu: \"%.120s%s\"",
                fp->current+1, fp->nparts, p->url,
                (strlen(p->url) > 120)? "..." : "");

            fp->currentfp = p->headers?
                  hopen(p->url, "r:", "httphdr:v", p->headers, NULL)
                : hopen(p->url, "r");

            if (fp->currentfp == NULL) return -1;
        }
        else return 0;  // No more parts, so we're truly at EOF
    }

    n = fp->currentfp->mobile?
          fp->currentfp->backend->read(fp->currentfp, buffer, nbytes)
        : hread(fp->currentfp, buffer, nbytes);

    if (n == 0) {
        // We're at EOF on this part, so set up the next part
        hFILE *prevfp = fp->currentfp;
        free_part(&fp->parts[fp->current]);
        fp->current++;
        fp->currentfp = NULL;
        if (hclose(prevfp) < 0) return -1;
        goto open_next;
    }

    return n;  // Number of bytes read by (or an error from) fp->currentfp
}

static ssize_t multipart_write(hFILE *fpv, const void *buffer, size_t nbytes)
{
    errno = EROFS;
    return -1;
}

static off_t multipart_seek(hFILE *fpv, off_t offset, int whence)
{
    errno = ESPIPE;
    return -1;
}

static int multipart_close(hFILE *fpv)
{
    hFILE_multipart *fp = (hFILE_multipart *) fpv;

    free_all_parts(fp);
    if (fp->currentfp) {
        if (hclose(fp->currentfp) < 0) return -1;
    }

    return 0;
}

static const struct hFILE_backend multipart_backend =
{
    multipart_read, multipart_write, multipart_seek, NULL, multipart_close
};

// Returns 'v' (valid value), 'i' (invalid; required GA4GH field missing),
// or upon encountering an unexpected token, that token's type.
// Explicit `return '?'` means a JSON parsing error, typically a member key
// that is not a string.  An unexpected token may be a valid token that was
// not the type expected for a particular GA4GH field, or it may be '?' or
// '\0' which should be propagated.
static char
parse_ga4gh_body_json(hFILE_multipart *fp, hFILE *json,
                      kstring_t *b, kstring_t *header)
{
    hts_json_token t;

    if (hts_json_fnext(json, &t, b) != '{') return t.type;
    while (hts_json_fnext(json, &t, b) != '}') {
        if (t.type != 's') return '?';

        if (strcmp(t.str, "urls") == 0) {
            if (hts_json_fnext(json, &t, b) != '[') return t.type;

            while (hts_json_fnext(json, &t, b) != ']') {
                hfile_part *part;
                size_t n = 0, max = 0;

                hts_expand(hfile_part, fp->nparts+1, fp->maxparts, fp->parts);
                part = &fp->parts[fp->nparts++];
                part->url = NULL;
                part->headers = NULL;

                if (t.type != '{') return t.type;
                while (hts_json_fnext(json, &t, b) != '}') {
                    if (t.type != 's') return '?';

                    if (strcmp(t.str, "url") == 0) {
                        if (hts_json_fnext(json, &t, b) != 's') return t.type;
                        part->url = ks_release(b);
                    }
                    else if (strcmp(t.str, "headers") == 0) {
                        if (hts_json_fnext(json, &t, b) != '{') return t.type;

                        while (hts_json_fnext(json, &t, header) != '}') {
                            if (t.type != 's') return '?';

                            if (hts_json_fnext(json, &t, b) != 's')
                                return t.type;

                            kputs(": ", header);
                            kputs(t.str, header);
                            n++;
                            hts_expand(char *, n+1, max, part->headers);
                            part->headers[n-1] = ks_release(header);
                            part->headers[n] = NULL;
                        }
                    }
                    else if (hts_json_fskip_value(json, '\0') != 'v')
                        return '?';
                }

                if (! part->url) return 'i';
            }
        }
        else if (strcmp(t.str, "format") == 0) {
            if (hts_json_fnext(json, &t, b) != 's') return t.type;

            hts_log_debug("GA4GH JSON redirection to multipart %s data", t.str);
        }
        else if (hts_json_fskip_value(json, '\0') != 'v') return '?';
    }

    return 'v';
}

// Returns 'v' (valid value), 'i' (invalid; required GA4GH field missing),
// or upon encountering an unexpected token, that token's type.
// Explicit `return '?'` means a JSON parsing error, typically a member key
// that is not a string.  An unexpected token may be a valid token that was
// not the type expected for a particular GA4GH field, or it may be '?' or
// '\0' which should be propagated.
static char
parse_ga4gh_redirect_json(hFILE_multipart *fp, hFILE *json,
                          kstring_t *b, kstring_t *header) {
    hts_json_token t;

    if (hts_json_fnext(json, &t, b) != '{') return t.type;
    while (hts_json_fnext(json, &t, b) != '}') {
        if (t.type != 's') return '?';

        if (strcmp(t.str, "htsget") == 0) {
            char ret = parse_ga4gh_body_json(fp, json, b, header);
            if (ret != 'v') return ret;
        }
        else return '?';
    }

    if (hts_json_fnext(json, &t, b) != '\0') return '?';

    return 'v';
}

hFILE *hopen_htsget_redirect(hFILE *hfile, const char *mode)
{
    hFILE_multipart *fp;
    kstring_t s1 = { 0, 0, NULL }, s2 = { 0, 0, NULL };
    char ret;

    fp = (hFILE_multipart *) hfile_init(sizeof (hFILE_multipart), mode, 0);
    if (fp == NULL) return NULL;

    fp->parts = NULL;
    fp->nparts = fp->maxparts = 0;

    ret = parse_ga4gh_redirect_json(fp, hfile, &s1, &s2);
    free(s1.s);
    free(s2.s);
    if (ret != 'v') {
        free_all_parts(fp);
        hfile_destroy((hFILE *) fp);
        errno = (ret == '?' || ret == '\0')? EPROTO : EINVAL;
        return NULL;
    }

    fp->current = 0;
    fp->currentfp = NULL;
    fp->base.backend = &multipart_backend;
    return &fp->base;
}
