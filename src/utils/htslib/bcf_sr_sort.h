/*
    Copyright (C) 2017 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

/*
    Reorder duplicate lines so that compatible variant types are
    returned together by bcf_sr_next_line()

    - readers grouped by variants. Even with many readers there will be
      typically only several groups

*/

#ifndef __BCF_SR_SORT_H__
#define __BCF_SR_SORT_H__

#include "htslib/synced_bcf_reader.h"
#include "htslib/kbitset.h"

typedef struct
{
    int nrec, mrec;
    bcf1_t **rec;
}
vcf_buf_t;

typedef struct
{
    char *str;      // "A>C" for biallelic records or "A>C,A>CC" for multiallelic records
    int type;       // VCF_SNP, VCF_REF, etc.
    int nalt;       // number of alternate alleles in this record
    int nvcf, mvcf, *vcf;   // the list of readers with the same variants
    bcf1_t **rec;           // list of VCF records in the readers
    kbitset_t *mask;        // which groups contain the variant
}
var_t;

typedef struct
{
    char *key;              // only for debugging
    int nvar, mvar, *var;   // the variants and their type
    int nvcf;               // number of readers with the same variants
}
grp_t;

typedef struct
{
    int nvar, mvar, *var;   // list of compatible variants that can be output together
    int cnt;                // number of readers in this group
    kbitset_t *mask;        // which groups are populated in this set (replace with expandable bitmask)
}
varset_t;

typedef struct
{
    uint8_t score[256];
    int nvar, mvar;
    var_t *var;             // list of all variants from all readers
    int nvset, mvset;
    int mpmat, *pmat;       // pairing matrix, i-th vset and j-th group accessible as i*ngrp+j
    int ngrp, mgrp;
    int mcnt, *cnt;         // number of VCF covered by a varset
    grp_t *grp;             // list of VCF representatives, each with a unique combination of duplicate lines
    varset_t *vset;         // list of variant sets - combinations of compatible variants across multiple groups ready for output
    vcf_buf_t *vcf_buf;     // records sorted in output order, for each VCF
    bcf_srs_t *sr;
    void *grp_str2int;
    void *var_str2int;
    kstring_t str;
    int moff, noff, *off, mcharp;
    char **charp;
    const char *chr;
    int pos, nsr, msr;
    int pair;
}
sr_sort_t;

sr_sort_t *bcf_sr_sort_init(sr_sort_t *srt);
int bcf_sr_sort_next(bcf_srs_t *readers, sr_sort_t *srt, const char *chr, int pos);
void bcf_sr_sort_destroy(sr_sort_t *srt);
void bcf_sr_sort_remove_reader(bcf_srs_t *readers, sr_sort_t *srt, int i);

#endif
