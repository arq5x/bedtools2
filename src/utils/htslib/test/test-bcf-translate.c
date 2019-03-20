/*  test/test-bcf-translate.c

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

#include <config.h>
#include <stdio.h>
#include <htslib/vcf.h>

int main(int argc, char **argv)
{
    char *fname = argc>1 ? argv[1] : "/dev/null";
    htsFile *fp = hts_open(fname, "w");
    bcf_hdr_t *hdr1, *hdr2;

    hdr1 = bcf_hdr_init("w");
    hdr2 = bcf_hdr_init("w");

    // Add two shared and two private annotations
    bcf_hdr_append(hdr1, "##contig=<ID=1>");
    bcf_hdr_append(hdr1, "##contig=<ID=2>");
    bcf_hdr_append(hdr2, "##contig=<ID=2>");
    bcf_hdr_append(hdr2, "##contig=<ID=1>");
    bcf_hdr_append(hdr1, "##FILTER=<ID=FLT1,Description=\"Filter 1\">");
    bcf_hdr_append(hdr1, "##FILTER=<ID=FLT2,Description=\"Filter 2\">");
    bcf_hdr_append(hdr1, "##FILTER=<ID=FLT3,Description=\"Filter 3\">");
    bcf_hdr_append(hdr2, "##FILTER=<ID=FLT4,Description=\"Filter 4\">");
    bcf_hdr_append(hdr2, "##FILTER=<ID=FLT3,Description=\"Filter 3\">");
    bcf_hdr_append(hdr2, "##FILTER=<ID=FLT2,Description=\"Filter 2\">");
    bcf_hdr_append(hdr1, "##INFO=<ID=INF1,Number=.,Type=Integer,Description=\"Info 1\">");
    bcf_hdr_append(hdr1, "##INFO=<ID=INF2,Number=.,Type=Integer,Description=\"Info 2\">");
    bcf_hdr_append(hdr1, "##INFO=<ID=INF3,Number=.,Type=Integer,Description=\"Info 3\">");
    bcf_hdr_append(hdr2, "##INFO=<ID=INF4,Number=.,Type=Integer,Description=\"Info 4\">");
    bcf_hdr_append(hdr2, "##INFO=<ID=INF3,Number=.,Type=Integer,Description=\"Info 3\">");
    bcf_hdr_append(hdr2, "##INFO=<ID=INF2,Number=.,Type=Integer,Description=\"Info 2\">");
    bcf_hdr_append(hdr1, "##FORMAT=<ID=FMT1,Number=.,Type=Integer,Description=\"FMT 1\">");
    bcf_hdr_append(hdr1, "##FORMAT=<ID=FMT2,Number=.,Type=Integer,Description=\"FMT 2\">");
    bcf_hdr_append(hdr1, "##FORMAT=<ID=FMT3,Number=.,Type=Integer,Description=\"FMT 3\">");
    bcf_hdr_append(hdr2, "##FORMAT=<ID=FMT4,Number=.,Type=Integer,Description=\"FMT 4\">");
    bcf_hdr_append(hdr2, "##FORMAT=<ID=FMT3,Number=.,Type=Integer,Description=\"FMT 3\">");
    bcf_hdr_append(hdr2, "##FORMAT=<ID=FMT2,Number=.,Type=Integer,Description=\"FMT 2\">");
    bcf_hdr_add_sample(hdr1,"SMPL1");
    bcf_hdr_add_sample(hdr1,"SMPL2");
    bcf_hdr_add_sample(hdr2,"SMPL1");
    bcf_hdr_add_sample(hdr2,"SMPL2");
    bcf_hdr_sync(hdr1);
    bcf_hdr_sync(hdr2);

    hdr2 = bcf_hdr_merge(hdr2,hdr1);
    bcf_hdr_sync(hdr2);
    bcf_hdr_write(fp, hdr2);

    bcf1_t *rec = bcf_init1();
    rec->rid = bcf_hdr_name2id(hdr1, "1");
    rec->pos = 0;
    bcf_update_alleles_str(hdr1, rec, "G,A");
    int32_t tmpi[3];
    tmpi[0] = bcf_hdr_id2int(hdr1, BCF_DT_ID, "FLT1");
    tmpi[1] = bcf_hdr_id2int(hdr1, BCF_DT_ID, "FLT2");
    tmpi[2] = bcf_hdr_id2int(hdr1, BCF_DT_ID, "FLT3");
    bcf_update_filter(hdr1, rec, tmpi, 3);
    tmpi[0] = 1; bcf_update_info_int32(hdr1, rec, "INF1", tmpi, 1);
    tmpi[0] = 2; bcf_update_info_int32(hdr1, rec, "INF2", tmpi, 1);
    tmpi[0] = 3; bcf_update_info_int32(hdr1, rec, "INF3", tmpi, 1);
    tmpi[0] = tmpi[1] = 1; bcf_update_format_int32(hdr1, rec, "FMT1", tmpi, 2);
    tmpi[0] = tmpi[1] = 2; bcf_update_format_int32(hdr1, rec, "FMT2", tmpi, 2);
    tmpi[0] = tmpi[1] = 3; bcf_update_format_int32(hdr1, rec, "FMT3", tmpi, 2);

    bcf_remove_filter(hdr1, rec, bcf_hdr_id2int(hdr1, BCF_DT_ID, "FLT2"), 0);
    bcf_update_info_int32(hdr1, rec, "INF2", NULL, 0);
    bcf_update_format_int32(hdr1, rec, "FMT2", NULL, 0);

    bcf_translate(hdr2, hdr1, rec);
    bcf_write(fp, hdr2, rec);

    // Clean
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr1);
    bcf_hdr_destroy(hdr2);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
    return 0;
}



    //      // Create VCF header
    //      kstring_t str = {0,0,0};
    //      bcf_hdr_add_sample(hdr, "NA00003");
    //      bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
    //      bcf_hdr_write(fp, hdr);
    //      // Add a record
    //      // 20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
    //      // .. CHROM
    //      rec->rid = bcf_hdr_name2id(hdr, "20");
    //      // .. POS
    //      rec->pos = 14369;
    //      // .. ID
    //      bcf_update_id(hdr, rec, "rs6054257");
    //      // .. REF and ALT
    //      bcf_update_alleles_str(hdr, rec, "G,A");
    //      // .. QUAL
    //      rec->qual = 29;
    //      // .. FILTER
    //      int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
    //      bcf_update_filter(hdr, rec, &tmpi, 1);
    //      // .. INFO
    //      tmpi = 3;
    //      bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
    //      tmpi = 14;
    //      bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
    //      float tmpf = 0.5;
    //      bcf_update_info_float(hdr, rec, "AF", &tmpf, 1);
    //      bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
    //      bcf_update_info_flag(hdr, rec, "H2", NULL, 1);
    //      // .. FORMAT
    //      int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
    //      tmpia[0] = bcf_gt_phased(0);
    //      tmpia[1] = bcf_gt_phased(0);
    //      tmpia[2] = bcf_gt_phased(1);
    //      tmpia[3] = bcf_gt_phased(0);
    //      tmpia[4] = bcf_gt_unphased(1);
    //      tmpia[5] = bcf_gt_unphased(1);
    //      bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
    //      tmpia[0] = 48;
    //      tmpia[1] = 48;
    //      tmpia[2] = 43;
    //      bcf_update_format_int32(hdr, rec, "GQ", tmpia, bcf_hdr_nsamples(hdr));
    //      tmpia[0] = 1;
    //      tmpia[1] = 8;
    //      tmpia[2] = 5;
    //      bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr));
    //      tmpia[0] = 51;
    //      tmpia[1] = 51;
    //      tmpia[2] = 51;
    //      tmpia[3] = 51;
    //      tmpia[4] = bcf_int32_missing;
    //      tmpia[5] = bcf_int32_missing;
    //      bcf_update_format_int32(hdr, rec, "HQ", tmpia, bcf_hdr_nsamples(hdr)*2);
    //      char *tmp_str[] = {"String1","SomeOtherString2","YetAnotherString3"};
    //      bcf_update_format_string(hdr, rec, "TS", (const char**)tmp_str, 3);
    //      bcf_write1(fp, hdr, rec);
    //      // 20     1110696 . A      G,T     67   .   NS=2;DP=10;AF=0.333,.;AA=T;DB GT 2 1   ./.
    //      bcf_clear1(rec);
    //      rec->rid = bcf_hdr_name2id(hdr, "20");
    //      rec->pos = 1110695;
    //      bcf_update_alleles_str(hdr, rec, "A,G,T");
    //      rec->qual = 67;
    //      tmpi = 2;
    //      bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
    //      tmpi = 10;
    //      bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
    //      float *tmpfa = (float*)malloc(2*sizeof(float));
    //      tmpfa[0] = 0.333;
    //      bcf_float_set_missing(tmpfa[1]);
    //      bcf_update_info_float(hdr, rec, "AF", tmpfa, 2);
    //      bcf_update_info_string(hdr, rec, "AA", "T");
    //      bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
    //      tmpia[0] = bcf_gt_phased(2);
    //      tmpia[1] = bcf_int32_vector_end;
    //      tmpia[2] = bcf_gt_phased(1);
    //      tmpia[3] = bcf_int32_vector_end;
    //      tmpia[4] = bcf_gt_missing;
    //      tmpia[5] = bcf_gt_missing;
    //      bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
    //      bcf_write1(fp, hdr, rec);
    //      free(tmpia);
    //      free(tmpfa);
