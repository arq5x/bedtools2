/* The MIT/Expat License

Copyright (C) 2017 Genome Research Ltd.

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
DEALINGS IN THE SOFTWARE.
 */
// Spam seeks
#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"

int main(int argc, char *argv[]) {
    if (argc <= 1) {
        fprintf(stderr, "Usage: thrash_threads4 input.bam\n");
        exit(1);
    }

    // Find a valid seek location ~64M into the file
    int i;
    ssize_t got;
    BGZF *fpin  = bgzf_open(argv[1], "r");
    uint64_t upos = 0, uend = 0;
    char buf[100000];
    for (i = 0; i < 100; i++) {
        if ((got = bgzf_read(fpin, buf, 65536)) < 0)
            abort();
        upos += got;
    }
    int64_t pos = bgzf_tell(fpin);
    while ((got = bgzf_read(fpin, buf, 65536)) > 0) {
        uend += got;
    }
    if (got < 0) abort();
    int64_t end = bgzf_tell(fpin);
    bgzf_close(fpin);

    // Ensure input is big enough to avoid case 3,4 below going off the end
    // of the file
    if (uend < upos + 10000000) {
        fprintf(stderr, "Please supply a bigger input file\n");
        exit(1);
    }

#define N 1000

    // Spam random seeks & reads
    for (i = 0; i < 1000; i++) {
        printf("i=%d\t", i);
        fpin  = bgzf_open(argv[1], "r");
        int j, eof = 0, mt = 0;
        for (j = 0; j < 80; j++) {
            int n = rand() % 7;
            putchar('0'+n); fflush(stdout);
            switch (n) {
            case 0: // start
                if (bgzf_seek(fpin, 0LL, SEEK_SET) < 0) puts("!");//abort();
                eof = 0;
                break;
            case 1: // mid
                if (bgzf_seek(fpin, pos, SEEK_SET) < 0) puts("!");//abort();
                eof = 0;
                break;
            case 2: // eof
                if (bgzf_seek(fpin, end, SEEK_SET) < 0) puts("!");//abort();
                eof = 1;
                break;
            case 3: case 4: {
                int l = rand()%(n==3?100000:100);
                if (bgzf_read(fpin, buf, l) != l*(1-eof)) abort();
                break;
            }
            case 5:
                usleep(N);
                break;
            case 6:
                if (!mt)
                    bgzf_mt(fpin, 8, 256);
                mt = 1;
                break;
            }
        }
        printf("\n");
        if (bgzf_close(fpin))
            abort();
    }

    return 0;
}
