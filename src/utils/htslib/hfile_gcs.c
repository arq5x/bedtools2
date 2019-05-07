/*  hfile_gcs.c -- Google Cloud Storage backend for low-level file streams.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "hfile_internal.h"
#ifdef ENABLE_PLUGINS
#include "version.h"
#endif

static hFILE *
gcs_rewrite(const char *gsurl, const char *mode, int mode_has_colon,
            va_list *argsp)
{
    const char *bucket, *path, *access_token;
    kstring_t mode_colon = { 0, 0, NULL };
    kstring_t url = { 0, 0, NULL };
    kstring_t auth_hdr = { 0, 0, NULL };
    hFILE *fp = NULL;

    // GCS URL format is gs[+SCHEME]://BUCKET/PATH

    if (gsurl[2] == '+') {
        bucket = strchr(gsurl, ':') + 1;
        kputsn(&gsurl[3], bucket - &gsurl[3], &url);
    }
    else {
        kputs("https:", &url);
        bucket = &gsurl[3];
    }
    while (*bucket == '/') kputc(*bucket++, &url);

    path = bucket + strcspn(bucket, "/?#");

    kputsn(bucket, path - bucket, &url);
    if (strchr(mode, 'r')) kputs(".storage-download", &url);
    else if (strchr(mode, 'w')) kputs(".storage-upload", &url);
    else kputs(".storage", &url);
    kputs(".googleapis.com", &url);

    kputs(path, &url);

    if (hts_verbose >= 8)
        fprintf(stderr, "[M::gcs_open] rewrote URL as %s\n", url.s);

    // TODO Find the access token in a more standard way
    access_token = getenv("GCS_OAUTH_TOKEN");

    if (access_token) {
        kputs("Authorization: Bearer ", &auth_hdr);
        kputs(access_token, &auth_hdr);
    }

    if (argsp || auth_hdr.l > 0 || mode_has_colon) {
        if (! mode_has_colon) {
            kputs(mode, &mode_colon);
            kputc(':', &mode_colon);
            mode = mode_colon.s;
        }

        fp = hopen(url.s, mode, "va_list", argsp,
                   "httphdr", (auth_hdr.l > 0)? auth_hdr.s : NULL, NULL);
    }
    else
        fp = hopen(url.s, mode);

    free(mode_colon.s);
    free(url.s);
    free(auth_hdr.s);
    return fp;
}

static hFILE *gcs_open(const char *url, const char *mode)
{
    return gcs_rewrite(url, mode, 0, NULL);
}

static hFILE *gcs_vopen(const char *url, const char *mode_colon, va_list args0)
{
    // Need to use va_copy() as we can only take the address of an actual
    // va_list object, not that of a parameter as its type may have decayed.
    va_list args;
    va_copy(args, args0);
    hFILE *fp = gcs_rewrite(url, mode_colon, 1, &args);
    va_end(args);
    return fp;
}

int PLUGIN_GLOBAL(hfile_plugin_init,_gcs)(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { gcs_open, hfile_always_remote, "Google Cloud Storage",
          2000 + 50, gcs_vopen
        };

#ifdef ENABLE_PLUGINS
    // Embed version string for examination via strings(1) or what(1)
    static const char id[] = "@(#)hfile_gcs plugin (htslib)\t" HTS_VERSION;
    if (hts_verbose >= 9)
        fprintf(stderr, "[M::hfile_gcs.init] version %s\n", strchr(id, '\t')+1);
#endif

    self->name = "Google Cloud Storage";
    hfile_add_scheme_handler("gs", &handler);
    hfile_add_scheme_handler("gs+http", &handler);
    hfile_add_scheme_handler("gs+https", &handler);
    return 0;
}
