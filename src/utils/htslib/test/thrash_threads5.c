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
// A basic 'zcat filename [N-threads]'

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"

#define N 1000
int main(int argc, char *argv[]) {
    char buf[N];
    ssize_t l, t = 0;

    if (argc < 2 || isatty(STDOUT_FILENO)) {
        fprintf(stderr,
                "Usage: thrash_threads5 input.bam num_threads | md5sum\n");
        exit(1);
    }

    BGZF *fpin  = bgzf_open(argv[1], "r");
    hts_tpool *p = NULL;
    if (argc > 2) {
        p = hts_tpool_init(atoi(argv[2]));
        bgzf_thread_pool(fpin,  p, 0);
    }
    int n = rand()%(N-1)+1;
    while ((l = bgzf_read(fpin, buf, n)) > 0) {
        if (l != write(STDOUT_FILENO, buf, l)) abort();
        t += l;
        if (l != n) {
            fprintf(stderr, "expected %d bytes, got %d\n", n, (int)l);
            break;
        }
        n = rand()%(N-1)+1;
    }
    fprintf(stderr, "close=%d\n", (int)bgzf_close(fpin));
    if (p) hts_tpool_destroy(p);

    fprintf(stderr, "wrote %d bytes\n", (int)t);

    return 0;
}
