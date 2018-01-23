/*  textutils.c -- non-bioinformatics utility routines for text etc.

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

#include "htslib/hfile.h"
#include "htslib/kstring.h"

#include "hts_internal.h"

static int dehex(char c)
{
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    else if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    else if (c >= '0' && c <= '9') return c - '0';
    else return -1;  // Hence dehex('\0') = -1
}

int hts_decode_percent(char *dest, size_t *destlen, const char *s)
{
    char *d = dest;
    int hi, lo;

    while (*s) {
        if (*s == '%' && (hi = dehex(s[1])) >= 0 && (lo = dehex(s[2])) >= 0) {
            *d++ = (hi << 4) | lo;
            s += 3;
        }
        else *d++ = *s++;
    }

    *d = '\0';
    *destlen = d - dest;
    return 0;
}

static int debase64(char c)
{
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    else if (c >= 'A' && c <= 'Z') return c - 'A';
    else if (c >= '0' && c <= '9') return c - '0' + 52;
    else if (c == '/') return 63;
    else if (c == '+') return 62;
    else return -1;  // Hence debase64('\0') = -1
}

size_t hts_base64_decoded_length(size_t len)
{
    size_t nquartets = (len + 2) / 4;
    return 3 * nquartets;
}

int hts_decode_base64(char *dest, size_t *destlen, const char *s)
{
    char *d = dest;
    int x0, x1, x2, x3;

    while (1) {
        x0 = debase64(*s++);
        x1 = (x0 >= 0)? debase64(*s++) : -1;
        x2 = (x1 >= 0)? debase64(*s++) : -1;
        x3 = (x2 >= 0)? debase64(*s++) : -1;
        if (x3 < 0) break;

        *d++ = (x0 << 2) | (x1 >> 4);
        *d++ = (x1 << 4) | (x2 >> 2);
        *d++ = (x2 << 6) | x3;
    }

    if (x1 >= 0) *d++ = (x0 << 2) | (x1 >> 4);
    if (x2 >= 0) *d++ = (x1 << 4) | (x2 >> 2);

    *destlen = d - dest;
    return 0;
}

static char *encode_utf8(char *s, unsigned x)
{
    if (x >= 0x10000) {
        *s++ = 0xF0 | (x >> 18);
        *s++ = 0x80 | ((x >> 12) & 0x3F);
        *s++ = 0x80 | ((x >> 6) & 0x3F);
        *s++ = 0x80 | (x & 0x3F);
    }
    else if (x >= 0x800) {
        *s++ = 0xE0 | (x >> 12);
        *s++ = 0x80 | ((x >> 6) & 0x3F);
        *s++ = 0x80 | (x & 0x3F);
    }
    else if (x >= 0x80) {
        *s++ = 0xC0 | (x >> 6);
        *s++ = 0x80 | (x & 0x3F);
    }
    else *s++ = x;

    return s;
}

static char *sscan_string(char *s)
{
    char *d = s;
    int d1, d2, d3, d4;

    for (;;) switch (*s) {
    case '\\':
        switch (s[1]) {
        case '\0': *d = '\0'; return s+1;
        case 'b': *d++ = '\b'; s += 2; break;
        case 'f': *d++ = '\f'; s += 2; break;
        case 'n': *d++ = '\n'; s += 2; break;
        case 'r': *d++ = '\r'; s += 2; break;
        case 't': *d++ = '\t'; s += 2; break;
        default:  *d++ = s[1]; s += 2; break;
        case 'u':
            if ((d1 = dehex(s[2])) >= 0 && (d2 = dehex(s[3])) >= 0 &&
                (d3 = dehex(s[4])) >= 0 && (d4 = dehex(s[5])) >= 0) {
                d = encode_utf8(d, d1 << 12 | d2 << 8 | d3 << 4 | d4);
                s += 6;
            }
            break;
        }
        break;

    case '"':
        *d = '\0';
        return s+1;

    case '\0':
        *d = '\0';
        return s;

    default:
        *d++ = *s++;
        break;
    }
}

static void fscan_string(hFILE *fp, kstring_t *d)
{
    int c, d1, d2, d3, d4;

    while ((c = hgetc(fp)) != EOF) switch (c) {
    case '\\':
        if ((c = hgetc(fp)) == EOF) return;
        switch (c) {
        case 'b': kputc('\b', d); break;
        case 'f': kputc('\f', d); break;
        case 'n': kputc('\n', d); break;
        case 'r': kputc('\r', d); break;
        case 't': kputc('\t', d); break;
        default:  kputc(c,    d); break;
        case 'u':
            if ((c = hgetc(fp)) != EOF && (d1 = dehex(c)) >= 0 &&
                (c = hgetc(fp)) != EOF && (d2 = dehex(c)) >= 0 &&
                (c = hgetc(fp)) != EOF && (d3 = dehex(c)) >= 0 &&
                (c = hgetc(fp)) != EOF && (d4 = dehex(c)) >= 0) {
                char buf[8];
                char *lim = encode_utf8(buf, d1 << 12 | d2 << 8 | d3 << 4 | d4);
                kputsn(buf, lim - buf, d);
            }
            break;
        }
        break;

    case '"':
        return;

    default:
        kputc(c, d);
        break;
    }
}

static char token_type(hts_json_token *token)
{
    const char *s = token->str;

    switch (*s) {
    case 'f':
        return (strcmp(s, "false") == 0)? 'b' : '?';
    case 'n':
        return (strcmp(s, "null") == 0)? '.' : '?';
    case 't':
        return (strcmp(s, "true") == 0)? 'b' : '?';
    case '-':
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
        return 'n';
    default:
        return '?';
    }
}

char hts_json_snext(char *str, size_t *state, hts_json_token *token)
{
    char *s = &str[*state >> 2];
    int hidden = *state & 3;

    if (hidden) {
        *state &= ~3;
        return token->type = "?}]?"[hidden];
    }

#define STATE(s,h)  (((s) - str) << 2 | (h))

    for (;;) switch (*s) {
    case ' ':
    case '\t':
    case '\r':
    case '\n':
    case ',':
    case ':':
        s++;
        continue;

    case '\0':
        return token->type = '\0';

    case '{':
    case '[':
    case '}':
    case ']':
        *state = STATE(s+1, 0);
        return token->type = *s;

    case '"':
        token->str = s+1;
        *state = STATE(sscan_string(s+1), 0);
        return token->type = 's';

    default:
        token->str = s;
        s += strcspn(s, " \t\r\n,]}");
        hidden = (*s == '}')? 1 : (*s == ']')? 2 : 0;
        if (*s != '\0') *s++ = '\0';
        *state = STATE(s, hidden);
        return token->type = token_type(token);
    }

#undef STATE
}

char hts_json_fnext(struct hFILE *fp, hts_json_token *token, kstring_t *kstr)
{
    char peek;
    int c;

    for (;;) switch (c = hgetc(fp)) {
    case ' ':
    case '\t':
    case '\r':
    case '\n':
    case ',':
    case ':':
        continue;

    case EOF:
        return token->type = '\0';

    case '{':
    case '[':
    case '}':
    case ']':
        return token->type = c;

    case '"':
        kstr->l = 0;
        fscan_string(fp, kstr);
        if (kstr->l == 0) kputsn("", 0, kstr);
        token->str = kstr->s;
        return token->type = 's';

    default:
        kstr->l = 0;
        kputc(c, kstr);
        while (hpeek(fp, &peek, 1) == 1 && !strchr(" \t\r\n,]}", peek)) {
            if ((c = hgetc(fp)) == EOF) break;
            kputc(c, kstr);
        }
        token->str = kstr->s;
        return token->type = token_type(token);
    }
}


typedef char hts_json_nextfn(void *arg1, void *arg2, hts_json_token *token);

static char skip_value(char type, hts_json_nextfn *next, void *arg1, void *arg2)
{
    hts_json_token token;
    int level;

    switch (type? type : next(arg1, arg2, &token)) {
    case '\0':
        return '\0';

    case '?':
    case '}':
    case ']':
        return '?';

    case '{':
    case '[':
        level = 1;
        break;

    default:
        return 'v';
    }

    while (level > 0)
        switch (next(arg1, arg2, &token)) {
        case '\0':
            return '\0';

        case '?':
            return '?';

        case '{':
        case '[':
            level++;
            break;

        case '}':
        case ']':
            --level;
            break;

        default:
            break;
        }

    return 'v';
}

static char snext(void *arg1, void *arg2, hts_json_token *token)
{
    return hts_json_snext(arg1, arg2, token);
}
char hts_json_sskip_value(char *str, size_t *state, char type)
{
    return skip_value(type, snext, str, state);
}

static char fnext(void *arg1, void *arg2, hts_json_token *token)
{
    return hts_json_fnext(arg1, token, arg2);
}
char hts_json_fskip_value(struct hFILE *fp, char type)
{
    kstring_t str = { 0, 0, NULL };
    char ret = skip_value(type, fnext, fp, &str);
    free(str.s);
    return ret;
}
