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
    BGZF *fpin  = bgzf_open(argv[1], "r");
    char buf[65536];
    for (i = 0; i < 1000; i++)
        if (bgzf_read(fpin, buf, 65536) < 0)
            abort();
    int64_t pos = bgzf_tell(fpin);
    bgzf_close(fpin);

#define N 1000

    // Spam seeks
    for (i = 0; i < 1000; i++) {
        printf("i=%d\n", i);
        fpin  = bgzf_open(argv[1], "r");
        bgzf_mt(fpin, 8, 256);
        if (bgzf_seek(fpin, pos, SEEK_SET) < 0) puts("!");//abort();
        usleep(N);
        //if (bgzf_read(fpin, buf, 65536) < 0) abort();
        //write(1, buf, 65536);
        if (bgzf_seek(fpin, 0LL, SEEK_SET) < 0) puts("!");//abort();
        usleep(N);
        //if (bgzf_read(fpin, buf, 65536) < 0) abort();
        //write(1, buf, 65536);
        if (bgzf_close(fpin))
            abort();
    }

    return 0;
}
