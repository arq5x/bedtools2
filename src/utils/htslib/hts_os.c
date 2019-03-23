/// @file hts_os.c
/// Operating System specific tweaks, for compatibility with POSIX.
/*
   Copyright (C) 2017 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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

// Windows (maybe more) lack a drand48 implementation.
#ifndef HAVE_DRAND48
#include "os/rand.c"
#else
#include <stdlib.h>
void hts_srand48(long seed) { srand48(seed); }
double hts_erand48(unsigned short xseed[3]) { return erand48(xseed); }
double hts_drand48(void) { return drand48(); }
double hts_lrand48(void) { return lrand48(); }
#endif

// // On Windows when using the MSYS or Cygwin terminals, isatty fails
// #ifdef _WIN32
// #define USE_FILEEXTD
// #include "os/iscygpty.c"
// #endif
