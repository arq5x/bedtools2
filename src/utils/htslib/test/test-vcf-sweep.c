/*  test/test-vcf-sweep.c -- VCF test harness.

    Copyright (C) 2013 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <stdio.h>
#include <htslib/vcf_sweep.h>

int main(int argc, char **argv)
{
    if ( argc!=2 )
    {
        fprintf(stderr,"Usage: test-vcf-sweep <file.bcf|file.vcf>\n");
        return 1;
    }

    // Init variables. The checksum is just for this test program to output
    // something and verify that all sites are read in both passes - fwd and
    // bwd.
    bcf_sweep_t *sw = bcf_sweep_init(argv[1]);
    bcf_hdr_t *hdr  = bcf_sweep_hdr(sw);
    int chksum = 0;

    // First we must sweep forward and read the whole file to build an index.
    // If this is undesirable, we can require the presence of a .gzi index
    // which can be created with `bgzip -r` from the samtools/htslib package
    bcf1_t *rec;
    while ( (rec = bcf_sweep_fwd(sw)) ) chksum += rec->pos+1;
    printf("fwd position chksum: %d\n", chksum);

    // Now sweep backward.
    chksum = 0;
    while ( (rec = bcf_sweep_bwd(sw)) ) chksum += rec->pos+1;
    printf("bwd position chksum: %d\n", chksum);

    // And forward and backward again, this time summing the PL vectors
    int i,j, mPLs = 0, nPLs;
    int32_t *PLs = NULL;
    chksum = 0;
    while ( (rec = bcf_sweep_fwd(sw)) )
    {
        // get copy of the PL vectors
        nPLs = bcf_get_format_int32(hdr, rec, "PL", &PLs, &mPLs);
        if ( !nPLs ) continue;  // PL not present

        // how many values are there per sample
        int nvals = nPLs / bcf_hdr_nsamples(hdr);

        int32_t *ptr = PLs;
        for (i=0; i<bcf_hdr_nsamples(hdr); i++)
        {
            for (j=0; j<nvals; j++)
            {
                // check for shorter vectors (haploid genotypes amongst diploids)
                if ( ptr[j]==bcf_int32_vector_end ) break;

                // skip missing values
                if ( ptr[j]==bcf_int32_missing ) continue;

                chksum += ptr[j];
            }
            ptr += nvals;
        }
    }
    printf("fwd PL chksum: %d\n", chksum);

    // And the same backwards..
    chksum = 0;
    while ( (rec = bcf_sweep_bwd(sw)) )
    {
        nPLs = bcf_get_format_int32(hdr, rec, "PL", &PLs, &mPLs);
        if ( !nPLs ) continue;
        int nvals = nPLs / bcf_hdr_nsamples(hdr);
        int32_t *ptr = PLs;
        for (i=0; i<bcf_hdr_nsamples(hdr); i++)
        {
            for (j=0; j<nvals; j++)
            {
                if ( ptr[j]==bcf_int32_vector_end ) break;
                if ( ptr[j]==bcf_int32_missing ) continue;
                chksum += ptr[j];
            }
            ptr += nvals;
        }
    }
    printf("bwd PL chksum: %d\n", chksum);

    // Clean up
    bcf_sweep_destroy(sw);
    return 0;
}


