/*  hfile_s3.c -- Amazon S3 backend for low-level file streams.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "hts_internal.h"
#include "hfile_internal.h"
#ifdef ENABLE_PLUGINS
#include "version.h"
#endif
#include "htslib/hts.h"  // for hts_version() and hts_verbose
#include "htslib/kstring.h"

typedef struct {
    kstring_t id;
    kstring_t token;
    kstring_t secret;
    char *bucket;
    kstring_t auth_hdr;
    time_t auth_time;
    char date[40];
    char mode;
    char *headers[3];
} s3_auth_data;

#define AUTH_LIFETIME 60

#if defined HAVE_COMMONCRYPTO

#include <CommonCrypto/CommonHMAC.h>

#define DIGEST_BUFSIZ CC_SHA1_DIGEST_LENGTH

static size_t
s3_sign(unsigned char *digest, kstring_t *key, kstring_t *message)
{
    CCHmac(kCCHmacAlgSHA1, key->s, key->l, message->s, message->l, digest);
    return CC_SHA1_DIGEST_LENGTH;
}

#elif defined HAVE_HMAC

#include <openssl/hmac.h>

#define DIGEST_BUFSIZ EVP_MAX_MD_SIZE

static size_t
s3_sign(unsigned char *digest, kstring_t *key, kstring_t *message)
{
    unsigned int len;
    HMAC(EVP_sha1(), key->s, key->l,
         (unsigned char *) message->s, message->l, digest, &len);
    return len;
}

#else
#error No HMAC() routine found by configure
#endif

static void
urldecode_kput(const char *s, int len, kstring_t *str)
{
    char buf[3];
    int i = 0;

    while (i < len)
        if (s[i] == '%' && i+2 < len) {
            buf[0] = s[i+1], buf[1] = s[i+2], buf[2] = '\0';
            kputc(strtol(buf, NULL, 16), str);
            i += 3;
        }
        else kputc(s[i++], str);
}

static void base64_kput(const unsigned char *data, size_t len, kstring_t *str)
{
    static const char base64[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    size_t i = 0;
    unsigned x = 0;
    int bits = 0, pad = 0;

    while (bits || i < len) {
        if (bits < 6) {
            x <<= 8, bits += 8;
            if (i < len) x |= data[i++];
            else pad++;
        }

        bits -= 6;
        kputc(base64[(x >> bits) & 63], str);
    }

    str->l -= pad;
    kputsn("==", pad, str);
}

static int is_dns_compliant(const char *s0, const char *slim)
{
    int has_nondigit = 0, len = 0;
    const char *s;

    for (s = s0; s < slim; len++, s++)
        if (islower_c(*s))
            has_nondigit = 1;
        else if (*s == '-') {
            has_nondigit = 1;
            if (s == s0 || s+1 == slim) return 0;
        }
        else if (isdigit_c(*s))
            ;
        else if (*s == '.') {
            if (s == s0 || ! isalnum_c(s[-1])) return 0;
            if (s+1 == slim || ! isalnum_c(s[1])) return 0;
        }
        else return 0;

    return has_nondigit && len >= 3 && len <= 63;
}

static FILE *expand_tilde_open(const char *fname, const char *mode)
{
    FILE *fp;

    if (strncmp(fname, "~/", 2) == 0) {
        kstring_t full_fname = { 0, 0, NULL };
        const char *home = getenv("HOME");
        if (! home) return NULL;

        kputs(home, &full_fname);
        kputs(&fname[1], &full_fname);

        fp = fopen(full_fname.s, mode);
        free(full_fname.s);
    }
    else
        fp = fopen(fname, mode);

    return fp;
}

static void parse_ini(const char *fname, const char *section, ...)
{
    kstring_t line = { 0, 0, NULL };
    int active = 1;  // Start active, so global properties are accepted
    char *s;

    FILE *fp = expand_tilde_open(fname, "r");
    if (fp == NULL) return;

    while (line.l = 0, kgetline(&line, (kgets_func *) fgets, fp) >= 0)
        if (line.s[0] == '[' && (s = strchr(line.s, ']')) != NULL) {
            *s = '\0';
            active = (strcmp(&line.s[1], section) == 0);
        }
        else if (active && (s = strpbrk(line.s, ":=")) != NULL) {
            const char *key = line.s, *value = &s[1], *akey;
            va_list args;

            while (isspace_c(*key)) key++;
            while (s > key && isspace_c(s[-1])) s--;
            *s = '\0';

            while (isspace_c(*value)) value++;
            while (line.l > 0 && isspace_c(line.s[line.l-1]))
                line.s[--line.l] = '\0';

            va_start(args, section);
            while ((akey = va_arg(args, const char *)) != NULL) {
                kstring_t *avar = va_arg(args, kstring_t *);
                if (strcmp(key, akey) == 0) { kputs(value, avar); break; }
            }
            va_end(args);
        }

    fclose(fp);
    free(line.s);
}

static void parse_simple(const char *fname, kstring_t *id, kstring_t *secret)
{
    kstring_t text = { 0, 0, NULL };
    char *s;
    size_t len;

    FILE *fp = expand_tilde_open(fname, "r");
    if (fp == NULL) return;

    while (kgetline(&text, (kgets_func *) fgets, fp) >= 0)
        kputc(' ', &text);
    fclose(fp);

    s = text.s;
    while (isspace_c(*s)) s++;
    kputsn(s, len = strcspn(s, " \t"), id);

    s += len;
    while (isspace_c(*s)) s++;
    kputsn(s, strcspn(s, " \t"), secret);

    free(text.s);
}

static int copy_auth_headers(s3_auth_data *ad, char ***hdrs) {
    char **hdr = &ad->headers[0];
    *hdrs = hdr;
    *hdr = strdup(ad->date);
    if (!*hdr) return -1;
    hdr++;
    if (ad->auth_hdr.l) {
        *hdr = strdup(ad->auth_hdr.s);
        if (!*hdr) { free(ad->headers[0]); return -1; }
        hdr++;
    }
    *hdr = NULL;
    return 0;
}

static void free_auth_data(s3_auth_data *ad) {
    free(ad->id.s);
    free(ad->token.s);
    free(ad->secret.s);
    free(ad->bucket);
    free(ad->auth_hdr.s);
    free(ad);
}

static int auth_header_callback(void *ctx, char ***hdrs) {
    s3_auth_data *ad = (s3_auth_data *) ctx;

    time_t now = time(NULL);
#ifdef HAVE_GMTIME_R
    struct tm tm_buffer;
    struct tm *tm = gmtime_r(&now, &tm_buffer);
#else
    struct tm *tm = gmtime(&now);
#endif
    kstring_t message = { 0, 0, NULL };
    unsigned char digest[DIGEST_BUFSIZ];
    size_t digest_len;

    if (!hdrs) { // Closing connection
        free_auth_data(ad);
        return 0;
    }

    if (now - ad->auth_time < AUTH_LIFETIME) {
        // Last auth string should still be valid
        *hdrs = NULL;
        return 0;
    }

    strftime(ad->date, sizeof(ad->date), "Date: %a, %d %b %Y %H:%M:%S GMT", tm);
    if (!ad->id.l || !ad->secret.l) {
        ad->auth_time = now;
        return copy_auth_headers(ad, hdrs);
    }

    if (ksprintf(&message, "%s\n\n\n%s\n%s%s%s/%s",
                 ad->mode == 'r' ? "GET" : "PUT", ad->date + 6,
                 ad->token.l ? "x-amz-security-token:" : "",
                 ad->token.l ? ad->token.s : "",
                 ad->token.l ? "\n" : "",
                 ad->bucket) < 0) {
        return -1;
    }

    digest_len = s3_sign(digest, &ad->secret, &message);
    ad->auth_hdr.l = 0;
    if (ksprintf(&ad->auth_hdr, "Authorization: AWS %s:", ad->id.s) < 0)
        goto fail;
    base64_kput(digest, digest_len, &ad->auth_hdr);

    free(message.s);
    ad->auth_time = now;
    return copy_auth_headers(ad, hdrs);

 fail:
    free(message.s);
    return -1;
}

static hFILE * s3_rewrite(const char *s3url, const char *mode, va_list *argsp)
{
    const char *bucket, *path;
    char *header_list[4], **header = header_list;

    kstring_t url = { 0, 0, NULL };
    kstring_t profile = { 0, 0, NULL };
    kstring_t host_base = { 0, 0, NULL };
    kstring_t token_hdr = { 0, 0, NULL };

    s3_auth_data *ad = calloc(1, sizeof(*ad));

    if (!ad)
        return NULL;
    ad->mode = strchr(mode, 'r') ? 'r' : 'w';

    // Our S3 URL format is s3[+SCHEME]://[ID[:SECRET[:TOKEN]]@]BUCKET/PATH

    if (s3url[2] == '+') {
        bucket = strchr(s3url, ':') + 1;
        kputsn(&s3url[3], bucket - &s3url[3], &url);
    }
    else {
        kputs("https:", &url);
        bucket = &s3url[3];
    }
    while (*bucket == '/') kputc(*bucket++, &url);

    path = bucket + strcspn(bucket, "/?#@");
    if (*path == '@') {
        const char *colon = strpbrk(bucket, ":@");
        if (*colon != ':') {
            urldecode_kput(bucket, colon - bucket, &profile);
        }
        else {
            const char *colon2 = strpbrk(&colon[1], ":@");
            urldecode_kput(bucket, colon - bucket, &ad->id);
            urldecode_kput(&colon[1], colon2 - &colon[1], &ad->secret);
            if (*colon2 == ':')
                urldecode_kput(&colon2[1], path - &colon2[1], &ad->token);
        }

        bucket = &path[1];
        path = bucket + strcspn(bucket, "/?#");
    }
    else {
        // If the URL has no ID[:SECRET]@, consider environment variables.
        const char *v;
        if ((v = getenv("AWS_ACCESS_KEY_ID")) != NULL) kputs(v, &ad->id);
        if ((v = getenv("AWS_SECRET_ACCESS_KEY")) != NULL) kputs(v, &ad->secret);
        if ((v = getenv("AWS_SESSION_TOKEN")) != NULL) kputs(v, &ad->token);

        if ((v = getenv("AWS_DEFAULT_PROFILE")) != NULL) kputs(v, &profile);
        else if ((v = getenv("AWS_PROFILE")) != NULL) kputs(v, &profile);
        else kputs("default", &profile);
    }

    if (ad->id.l == 0) {
        const char *v = getenv("AWS_SHARED_CREDENTIALS_FILE");
        parse_ini(v? v : "~/.aws/credentials", profile.s,
                  "aws_access_key_id", &ad->id,
                  "aws_secret_access_key", &ad->secret,
                  "aws_session_token", &ad->token, NULL);
    }
    if (ad->id.l == 0)
        parse_ini("~/.s3cfg", profile.s, "access_key", &ad->id,
                  "secret_key", &ad->secret, "access_token", &ad->token,
                  "host_base", &host_base, NULL);
    if (ad->id.l == 0)
        parse_simple("~/.awssecret", &ad->id, &ad->secret);

    if (host_base.l == 0)
        kputs("s3.amazonaws.com", &host_base);
    // Use virtual hosted-style access if possible, otherwise path-style.
    if (is_dns_compliant(bucket, path)) {
        kputsn(bucket, path - bucket, &url);
        kputc('.', &url);
        kputs(host_base.s, &url);
    }
    else {
        kputs(host_base.s, &url);
        kputc('/', &url);
        kputsn(bucket, path - bucket, &url);
    }
    kputs(path, &url);

    if (ad->token.l > 0) {
        kputs("X-Amz-Security-Token: ", &token_hdr);
        kputs(ad->token.s, &token_hdr);
        *header++ = token_hdr.s;
    }

    ad->bucket = strdup(bucket);
    if (!ad->bucket)
        goto fail;

    *header = NULL;
    hFILE *fp = hopen(url.s, mode, "va_list", argsp, "httphdr:v", header_list,
                      "httphdr_callback", auth_header_callback,
                      "httphdr_callback_data", ad, NULL);
    if (!fp) goto fail;

    free(url.s);
    free(profile.s);
    free(host_base.s);
    free(token_hdr.s);
    return fp;

 fail:
    free(url.s);
    free(profile.s);
    free(host_base.s);
    free(token_hdr.s);
    free_auth_data(ad);
    return NULL;
}

static hFILE *s3_open(const char *url, const char *mode)
{
    kstring_t mode_colon = { 0, 0, NULL };
    kputs(mode, &mode_colon);
    kputc(':', &mode_colon);
    hFILE *fp = s3_rewrite(url, mode_colon.s, NULL);
    free(mode_colon.s);
    return fp;
}

static hFILE *s3_vopen(const char *url, const char *mode_colon, va_list args0)
{
    // Need to use va_copy() as we can only take the address of an actual
    // va_list object, not that of a parameter whose type may have decayed.
    va_list args;
    va_copy(args, args0);
    hFILE *fp = s3_rewrite(url, mode_colon, &args);
    va_end(args);
    return fp;
}

int PLUGIN_GLOBAL(hfile_plugin_init,_s3)(struct hFILE_plugin *self)
{
    static const struct hFILE_scheme_handler handler =
        { s3_open, hfile_always_remote, "Amazon S3", 2000 + 50, s3_vopen
        };

#ifdef ENABLE_PLUGINS
    // Embed version string for examination via strings(1) or what(1)
    static const char id[] = "@(#)hfile_s3 plugin (htslib)\t" HTS_VERSION;
    if (hts_verbose >= 9)
        fprintf(stderr, "[M::hfile_s3.init] version %s\n", strchr(id, '\t')+1);
#endif

    self->name = "Amazon S3";
    hfile_add_scheme_handler("s3", &handler);
    hfile_add_scheme_handler("s3+http", &handler);
    hfile_add_scheme_handler("s3+https", &handler);
    return 0;
}
