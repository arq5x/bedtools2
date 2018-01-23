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
    Test bcf synced reader allele pairing
*/

#include <config.h>

#include <getopt.h>
#include <stdio.h>
#include <stdarg.h>
#include <htslib/synced_bcf_reader.h>

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

void usage(void)
{
    fprintf(stderr, "Usage: test-bcf-sr [OPTIONS] vcf-list.txt\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -p, --pair <logic[+ref]>     logic: snps,indels,both,snps+ref,indels+ref,both+ref,exact,some,all\n");
    fprintf(stderr, "\n");
    exit(-1);
}

int main(int argc, char *argv[])
{
    static struct option loptions[] =
    {
        {"help",no_argument,NULL,'h'},
        {"pair",required_argument,NULL,'p'},
        {NULL,0,NULL,0}
    };

    int c, pair = 0;
    while ((c = getopt_long(argc, argv, "p:h", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'p':
                if ( !strcmp(optarg,"snps") )            pair |= BCF_SR_PAIR_SNPS;
                else if ( !strcmp(optarg,"snp+ref") )    pair |= BCF_SR_PAIR_SNPS|BCF_SR_PAIR_SNP_REF;
                else if ( !strcmp(optarg,"snps+ref") )   pair |= BCF_SR_PAIR_SNPS|BCF_SR_PAIR_SNP_REF;
                else if ( !strcmp(optarg,"indels") )     pair |= BCF_SR_PAIR_INDELS;
                else if ( !strcmp(optarg,"indel+ref") )  pair |= BCF_SR_PAIR_INDELS|BCF_SR_PAIR_INDEL_REF;
                else if ( !strcmp(optarg,"indels+ref") ) pair |= BCF_SR_PAIR_INDELS|BCF_SR_PAIR_INDEL_REF;
                else if ( !strcmp(optarg,"both") )       pair |= BCF_SR_PAIR_BOTH;
                else if ( !strcmp(optarg,"both+ref") )   pair |= BCF_SR_PAIR_BOTH_REF;
                else if ( !strcmp(optarg,"any") )        pair |= BCF_SR_PAIR_ANY;
                else if ( !strcmp(optarg,"all") )        pair |= BCF_SR_PAIR_ANY;
                else if ( !strcmp(optarg,"some") )       pair |= BCF_SR_PAIR_SOME;
                else if ( !strcmp(optarg,"exact") )      pair  = BCF_SR_PAIR_EXACT;
                else error("The --pair logic \"%s\" not recognised.\n", optarg);
                break;
            default: usage();
        }
    }
    if ( !pair ) pair = BCF_SR_PAIR_EXACT;
    if ( optind == argc ) usage();

    int i, j, n, nvcf;
    char **vcf = hts_readlist(argv[optind], 1, &nvcf);
    if ( !vcf ) error("Could not parse %s\n", argv[optind]);

    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, pair);
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    for (i=0; i<nvcf; i++)
        if ( !bcf_sr_add_reader(sr,vcf[i]) ) error("Failed to open %s: %s\n", vcf[i],bcf_sr_strerror(sr->errnum));

    kstring_t str = {0,0,0};
    while ( (n=bcf_sr_next_line(sr)) )
    {
        for (i=0; i<sr->nreaders; i++)
        {
            if ( !bcf_sr_has_line(sr,i) ) continue;
            bcf1_t *rec = bcf_sr_get_line(sr, i);
            printf("%s:%d", bcf_seqname(bcf_sr_get_header(sr,i),rec),rec->pos+1);
            break;
        }

        for (i=0; i<sr->nreaders; i++)
        {
            printf("\t");

            if ( !bcf_sr_has_line(sr,i) )
            {
                printf("%s","-");
                continue;
            }

            str.l = 0;
            bcf1_t *rec = bcf_sr_get_line(sr, i);
            kputs(rec->n_allele > 1 ? rec->d.allele[1] : ".", &str);
            for (j=2; j<rec->n_allele; j++)
            {
                kputc(',', &str);
                kputs(rec->d.allele[j], &str);
            }
            printf("%s",str.s);
        }
        printf("\n");
    }

    free(str.s);
    bcf_sr_destroy(sr);
    for (i=0; i<nvcf; i++)
        free(vcf[i]);
    free(vcf);

    return 0;
}

