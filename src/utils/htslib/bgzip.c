/* bgzip.c -- Block compression/decompression utility.

   Copyright (C) 2008, 2009 Broad Institute / Massachusetts Institute of Technology
   Copyright (C) 2010, 2013-2018 Genome Research Ltd.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notices and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <getopt.h>
#include <inttypes.h>
#include <sys/stat.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

static const int WINDOW_SIZE = 64 * 1024;

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}

static int confirm_overwrite(const char *fn)
{
    int save_errno = errno;
    int ret = 0;

    if (isatty(STDIN_FILENO)) {
        char c;
        fprintf(stderr, "[bgzip] %s already exists; do you wish to overwrite (y or n)? ", fn);
        if (scanf("%c", &c) == 1 && (c == 'Y' || c == 'y')) ret = 1;
    }

    errno = save_errno;
    return ret;
}

static int bgzip_main_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Version: %s\n", hts_version());
    fprintf(stderr, "Usage:   bgzip [OPTIONS] [FILE] ...\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -b, --offset INT           decompress at virtual file pointer (0-based uncompressed offset)\n");
    fprintf(stderr, "   -c, --stdout               write on standard output, keep original files unchanged\n");
    fprintf(stderr, "   -d, --decompress           decompress\n");
    fprintf(stderr, "   -f, --force                overwrite files without asking\n");
    fprintf(stderr, "   -h, --help                 give this help\n");
    fprintf(stderr, "   -i, --index                compress and create BGZF index\n");
    fprintf(stderr, "   -I, --index-name FILE      name of BGZF index file [file.gz.gzi]\n");
    fprintf(stderr, "   -l, --compress-level INT   Compression level to use when compressing; 0 to 9, or -1 for default [-1]\n");
    fprintf(stderr, "   -r, --reindex              (re)index compressed file\n");
    fprintf(stderr, "   -g, --rebgzip              use an index file to bgzip a file\n");
    fprintf(stderr, "   -s, --size INT             decompress INT bytes (uncompressed size)\n");
    fprintf(stderr, "   -@, --threads INT          number of compression threads to use [1]\n");
    fprintf(stderr, "   -t, --test                 test integrity of compressed file");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char **argv)
{
    int c, compress, compress_level = -1, pstdout, is_forced, test, index = 0, rebgzip = 0, reindex = 0;
    BGZF *fp;
    void *buffer;
    long start, end, size;
    char *index_fname = NULL;
    int threads = 1;

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"offset", required_argument, NULL, 'b'},
        {"stdout", no_argument, NULL, 'c'},
        {"decompress", no_argument, NULL, 'd'},
        {"force", no_argument, NULL, 'f'},
        {"index", no_argument, NULL, 'i'},
        {"index-name", required_argument, NULL, 'I'},
        {"compress-level", required_argument, NULL, 'l'},
        {"reindex", no_argument, NULL, 'r'},
        {"rebgzip",no_argument,NULL,'g'},
        {"size", required_argument, NULL, 's'},
        {"threads", required_argument, NULL, '@'},
        {"test", no_argument, NULL, 't'},
        {"version", no_argument, NULL, 1},
        {NULL, 0, NULL, 0}
    };

    compress = 1; pstdout = 0; start = 0; size = -1; end = -1; is_forced = 0; test = 0;
    while((c  = getopt_long(argc, argv, "cdh?fb:@:s:iI:l:grt",loptions,NULL)) >= 0){
        switch(c){
        case 'd': compress = 0; break;
        case 'c': pstdout = 1; break;
        case 'b': start = atol(optarg); compress = 0; pstdout = 1; break;
        case 's': size = atol(optarg); pstdout = 1; break;
        case 'f': is_forced = 1; break;
        case 'i': index = 1; break;
        case 'I': index_fname = optarg; break;
        case 'l': compress_level = atol(optarg); break;
        case 'g': rebgzip = 1; break;
        case 'r': reindex = 1; compress = 0; break;
        case '@': threads = atoi(optarg); break;
        case 't': test = 1; compress = 0; reindex = 0; break;
        case 1:
            printf(
"bgzip (htslib) %s\n"
"Copyright (C) 2018 Genome Research Ltd.\n", hts_version());
            return EXIT_SUCCESS;
        case 'h':
        case '?': return bgzip_main_usage();
        }
    }
    if (size >= 0) end = start + size;
    if (end >= 0 && end < start) {
        fprintf(stderr, "[bgzip] Illegal region: [%ld, %ld]\n", start, end);
        return 1;
    }
    if (compress == 1) {
        struct stat sbuf;
        int f_src = fileno(stdin);
        char out_mode[3] = "w\0";
        char out_mode_exclusive[4] = "wx\0";

        if (compress_level < -1 || compress_level > 9) {
            fprintf(stderr, "[bgzip] Invalid compress-level: %d\n", compress_level);
            return 1;
        }
        if (compress_level >= 0) {
            out_mode[1] = compress_level + '0';
            out_mode_exclusive[2] = compress_level + '0';
        }

        if ( argc>optind )
        {
            if ( stat(argv[optind],&sbuf)<0 )
            {
                fprintf(stderr, "[bgzip] %s: %s\n", strerror(errno), argv[optind]);
                return 1;
            }

            if ((f_src = open(argv[optind], O_RDONLY)) < 0) {
                fprintf(stderr, "[bgzip] %s: %s\n", strerror(errno), argv[optind]);
                return 1;
            }

            if (pstdout)
                fp = bgzf_open("-", out_mode);
            else
            {
                char *name = malloc(strlen(argv[optind]) + 5);
                strcpy(name, argv[optind]);
                strcat(name, ".gz");
                fp = bgzf_open(name, is_forced? out_mode : out_mode_exclusive);
                if (fp == NULL && errno == EEXIST && confirm_overwrite(name))
                    fp = bgzf_open(name, out_mode);
                if (fp == NULL) {
                    fprintf(stderr, "[bgzip] can't create %s: %s\n", name, strerror(errno));
                    free(name);
                    return 1;
                }
                free(name);
            }
        }
        else if (!pstdout && isatty(fileno((FILE *)stdout)) )
            return bgzip_main_usage();
        else if ( index && !index_fname )
        {
            fprintf(stderr, "[bgzip] Index file name expected when writing to stdout\n");
            return 1;
        }
        else
            fp = bgzf_open("-", out_mode);

        if ( index && rebgzip )
        {
            fprintf(stderr, "[bgzip] Can't produce a index and rebgzip simultaneously\n");
            return 1;
        }

        if ( rebgzip && !index_fname )
        {
            fprintf(stderr, "[bgzip] Index file name expected when writing to stdout\n");
            return 1;
        }

        if (threads > 1)
            bgzf_mt(fp, threads, 256);

        if ( index ) bgzf_index_build_init(fp);
        buffer = malloc(WINDOW_SIZE);
#ifdef _WIN32
        _setmode(f_src, O_BINARY);
#endif
        if (rebgzip){
            if ( bgzf_index_load(fp, index_fname, NULL) < 0 ) error("Could not load index: %s.gzi\n", argv[optind]);

            while ((c = read(f_src, buffer, WINDOW_SIZE)) > 0)
                if (bgzf_block_write(fp, buffer, c) < 0) error("Could not write %d bytes: Error %d\n", c, fp->errcode);
        }
        else {
            while ((c = read(f_src, buffer, WINDOW_SIZE)) > 0)
                if (bgzf_write(fp, buffer, c) < 0) error("Could not write %d bytes: Error %d\n", c, fp->errcode);
        }
        if ( index )
        {
            if (index_fname) {
                if (bgzf_index_dump(fp, index_fname, NULL) < 0)
                    error("Could not write index to '%s'\n", index_fname);
            } else {
                if (bgzf_index_dump(fp, argv[optind], ".gz.gzi") < 0)
                    error("Could not write index to '%s.gz.gzi'", argv[optind]);
            }
        }
        if (bgzf_close(fp) < 0) error("Close failed: Error %d", fp->errcode);
        if (argc > optind && !pstdout) unlink(argv[optind]);
        free(buffer);
        close(f_src);
        return 0;
    }
    else if ( reindex )
    {
        if ( argc>optind )
        {
            fp = bgzf_open(argv[optind], "r");
            if ( !fp ) error("[bgzip] Could not open file: %s\n", argv[optind]);
        }
        else
        {
            if ( !index_fname ) error("[bgzip] Index file name expected when reading from stdin\n");
            fp = bgzf_open("-", "r");
            if ( !fp ) error("[bgzip] Could not read from stdin: %s\n", strerror(errno));
        }

        buffer = malloc(BGZF_BLOCK_SIZE);
        bgzf_index_build_init(fp);
        int ret;
        while ( (ret=bgzf_read(fp, buffer, BGZF_BLOCK_SIZE))>0 ) ;
        free(buffer);
        if ( ret<0 ) error("Is the file gzipped or bgzipped? The latter is required for indexing.\n");

        if ( index_fname ) {
            if (bgzf_index_dump(fp, index_fname, NULL) < 0)
                error("Could not write index to '%s'\n", index_fname);
        } else {
            if (bgzf_index_dump(fp, argv[optind], ".gzi") < 0)
                error("Could not write index to '%s.gzi'\n", argv[optind]);
        }

        if ( bgzf_close(fp)<0 ) error("Close failed: Error %d\n",fp->errcode);
        return 0;
    }
    else
    {
        struct stat sbuf;
        int f_dst;

        if ( argc>optind )
        {
            if ( stat(argv[optind],&sbuf)<0 )
            {
                fprintf(stderr, "[bgzip] %s: %s\n", strerror(errno), argv[optind]);
                return 1;
            }
            char *name;
            int len = strlen(argv[optind]);
            if ( strcmp(argv[optind]+len-3,".gz") && !test)
            {
                fprintf(stderr, "[bgzip] %s: unknown suffix -- ignored\n", argv[optind]);
                return 1;
            }
            fp = bgzf_open(argv[optind], "r");
            if (fp == NULL) {
                fprintf(stderr, "[bgzip] Could not open file: %s\n", argv[optind]);
                return 1;
            }

            if (pstdout || test) {
                f_dst = fileno(stdout);
            }
            else {
                const int wrflags = O_WRONLY | O_CREAT | O_TRUNC;
                name = strdup(argv[optind]);
                name[strlen(name) - 3] = '\0';
                f_dst = open(name, is_forced? wrflags : wrflags|O_EXCL, 0666);
                if (f_dst < 0 && errno == EEXIST && confirm_overwrite(name))
                    f_dst = open(name, wrflags, 0666);
                if (f_dst < 0) {
                    fprintf(stderr, "[bgzip] can't create %s: %s\n", name, strerror(errno));
                    free(name);
                    return 1;
                }
                free(name);
            }
        }
        else if (!pstdout && isatty(fileno((FILE *)stdin)) )
            return bgzip_main_usage();
        else
        {
            f_dst = fileno(stdout);
            fp = bgzf_open("-", "r");
            if (fp == NULL) {
                fprintf(stderr, "[bgzip] Could not read from stdin: %s\n", strerror(errno));
                return 1;
            }
        }

        if (!fp->is_compressed) {
            fprintf(stderr, "[bgzip] Expected compressed file -- ignored\n");
            return 1;
        }

        if (threads > 1)
            bgzf_mt(fp, threads, 256);

        buffer = malloc(WINDOW_SIZE);
        if ( start>0 )
        {
            if ( bgzf_index_load(fp, argv[optind], ".gzi") < 0 ) error("Could not load index: %s.gzi\n", argv[optind]);
            if ( bgzf_useek(fp, start, SEEK_SET) < 0 ) error("Could not seek to %d-th (uncompressd) byte\n", start);
        }
#ifdef _WIN32
        _setmode(f_dst, O_BINARY);
#endif
        while (1) {
            if (end < 0) c = bgzf_read(fp, buffer, WINDOW_SIZE);
            else c = bgzf_read(fp, buffer, (end - start > WINDOW_SIZE)? WINDOW_SIZE:(end - start));
            if (c == 0) break;
            if (c < 0) error("Error %d in block starting at offset %" PRId64 "(%" PRIX64 ")\n", fp->errcode, fp->block_address, fp->block_address);
            start += c;
            if ( !test && write(f_dst, buffer, c) != c ) {
#ifdef _WIN32
                if (GetLastError() != ERROR_NO_DATA)
#endif
                error("Could not write %d bytes\n", c);
            }
            if (end >= 0 && start >= end) break;
        }
        free(buffer);
        if (bgzf_close(fp) < 0) error("Close failed: Error %d\n",fp->errcode);
        if (!pstdout && !test) unlink(argv[optind]);
        return 0;
    }
}
