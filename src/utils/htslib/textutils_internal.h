/* textutils_internal.h -- non-bioinformatics utility routines for text etc.

   Copyright (C) 2016,2018 Genome Research Ltd.

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

#ifndef HTSLIB_TEXTUTILS_INTERNAL_H
#define HTSLIB_TEXTUTILS_INTERNAL_H

/* N.B. These interfaces may be used by plug-ins */

#include <ctype.h>
#include <stdlib.h>
#include "htslib/kstring.h"

#ifdef __cplusplus
extern "C" {
#endif

/// Decode percent-encoded (URL-encoded) text
/** On input, _dest_ should be a buffer at least the same size as _s_,
    and may be equal to _s_ to decode in place.  On output, _dest_ will be
    NUL-terminated and the number of characters written (not including the
    NUL) is stored in _destlen_.
*/
int hts_decode_percent(char *dest, size_t *destlen, const char *s);

/// Return decoded data length given length of base64-encoded text
/** This gives an upper bound, as it overestimates by a byte or two when
    the encoded text ends with (possibly omitted) `=` padding characters.
*/
size_t hts_base64_decoded_length(size_t len);

/// Decode base64-encoded data
/** On input, _dest_ should be a sufficient buffer (see `hts_base64_length()`),
    and may be equal to _s_ to decode in place.  On output, the number of
    bytes writen is stored in _destlen_.
*/
int hts_decode_base64(char *dest, size_t *destlen, const char *s);

/// Token structure returned by JSON lexing functions
/** Structure is defined in hts_internal.h
 */

typedef struct hts_json_token hts_json_token;

/// Allocate an empty JSON token structure, for use with hts_json_* functions
/** @return An empty token on success; NULL on failure
 */
hts_json_token * hts_json_alloc_token();

/// Free a JSON token
void hts_json_free_token(hts_json_token *token);

/// Accessor funtion to get JSON token type
/** @param  token Pointer to JSON token
    @return Character indicating the token type

Token types correspond to scalar JSON values and selected punctuation
as follows:
  - `s` string
  - `n` number
  - `b` boolean literal
  - `.` null literal
  - `{`, `}`, `[`, `]` object and array delimiters
  - `?` lexing error
  - `!` other errors (e.g. out of memory)
  - `\0` terminator at end of input
*/
char hts_json_token_type(hts_json_token *token);

/// Accessor funtion to get JSON token in string form
/** @param  token Pointer to JSON token
    @return String representation of the JSON token; NULL if unset

If the token was parsed from a string using hts_json_snext(), the return value
will point into the string passed as the first parameter to hts_json_snext().
If the token was parsed from a file using hts_json_fnext(), the return value
will point at the kstring_t buffer passed as the third parameter to
hts_json_fnext().  In that case, the value will only be valid until the
next call to hts_json_fnext().
 */
char *hts_json_token_str(hts_json_token *token);

/// Read one JSON token from a string
/** @param str    The input C string
    @param state  The input string state
    @param token  On return, filled in with the token read
    @return  The type of the token read

On return, `token->str` points into the supplied input string, which
is modified by having token-terminating characters overwritten as NULs.
The `state` argument records the current position within `str` after each
`hts_json_snext()` call, and should be set to 0 before the first call.
*/
char hts_json_snext(char *str, size_t *state, hts_json_token *token);

/// Read and discard a complete JSON value from a string
/** @param str    The input C string
    @param state  The input string state, as per `hts_json_snext()`
    @param type   If the first token of the value to be discarded has already
                  been read, provide its type; otherwise `'\0'`
    @return  One of `v` (success), `\0` (end of string), and `?` (lexing error)

Skips a complete JSON value, which may be a single token or an entire object
or array.
*/
char hts_json_sskip_value(char *str, size_t *state, char type);

/// Read one JSON token from a file
/** @param fp     The file stream
    @param token  On return, filled in with the token read
    @param kstr   Buffer used to store the token string returned
    @return  The type of the token read

The `kstr` buffer is used to store the string value of the token read,
so `token->str` is only valid until the next time `hts_json_fnext()` is
called with the same `kstr` argument.
*/
char hts_json_fnext(struct hFILE *fp, hts_json_token *token, kstring_t *kstr);

/// Read and discard a complete JSON value from a file
/** @param fp    The file stream
    @param type  If the first token of the value to be discarded has already
                 been read, provide its type; otherwise `'\0'`
    @return  One of `v` (success), `\0` (EOF), and `?` (lexing error)

Skips a complete JSON value, which may be a single token or an entire object
or array.
*/
char hts_json_fskip_value(struct hFILE *fp, char type);


// The <ctype.h> functions operate on ints such as are returned by fgetc(),
// i.e., characters represented as unsigned-char-valued ints, or EOF.
// To operate on plain chars (and to avoid warnings on some platforms),
// technically one must cast to unsigned char everywhere (see CERT STR37-C)
// or less painfully use these *_c() functions that operate on plain chars
// (but not EOF, which must be considered separately where it is applicable).
// TODO We may eventually wish to implement these functions directly without
// using their <ctype.h> equivalents, and thus make them immune to locales.
static inline int isalnum_c(char c) { return isalnum((unsigned char) c); }
static inline int isalpha_c(char c) { return isalpha((unsigned char) c); }
static inline int isdigit_c(char c) { return isdigit((unsigned char) c); }
static inline int isgraph_c(char c) { return isgraph((unsigned char) c); }
static inline int islower_c(char c) { return islower((unsigned char) c); }
static inline int isprint_c(char c) { return isprint((unsigned char) c); }
static inline int isspace_c(char c) { return isspace((unsigned char) c); }
static inline int isupper_c(char c) { return isupper((unsigned char) c); }
static inline char tolower_c(char c) { return tolower((unsigned char) c); }
static inline char toupper_c(char c) { return toupper((unsigned char) c); }

#ifdef __cplusplus
}
#endif

#endif
