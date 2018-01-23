/*  rand.h -- drand48 implementation from the FreeBSD source tree. */

/*
 * Copyright (c) 1993 Martin Birgmeier
 * All rights reserved.
 *
 * You may redistribute unmodified or modified versions of this source
 * code provided that the above copyright notice and this and the
 * following conditions are retained.
 *
 * This software is provided ``as is'', and comes with no warranties
 * of any kind. I shall in no event be liable for anything that happens
 * to anyone/anything when using this software.
 */

#ifndef HTSLIB_HTS_RAND_H
#define HTSLIB_HTS_RAND_H

void hts_srand48(long seed);
double hts_erand48(unsigned short xseed[3]);
double hts_drand48(void);
long hts_lrand48(void);

#endif /* HTSLIB_HTS_RAND_H */
