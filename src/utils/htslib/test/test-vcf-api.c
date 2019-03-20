/*  test/test-vcf-api.c -- VCF test harness.

    Copyright (C) 2013, 2014 Genome Research Ltd.

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
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

void write_bcf(char *fname)
{
    // Init
    htsFile *fp    = hts_open(fname,"wb");
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    bcf1_t *rec    = bcf_init1();

    // Create VCF header
    kstring_t str = {0,0,0};
    bcf_hdr_append(hdr, "##fileDate=20090805");
    bcf_hdr_append(hdr, "##FORMAT=<ID=UF,Number=1,Type=Integer,Description=\"Unused FORMAT\">");
    bcf_hdr_append(hdr, "##INFO=<ID=UI,Number=1,Type=Integer,Description=\"Unused INFO\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=Flt,Description=\"Unused FILTER\">");
    bcf_hdr_append(hdr, "##unused=<XX=AA,Description=\"Unused generic\">");
    bcf_hdr_append(hdr, "##unused=unformatted text 1");
    bcf_hdr_append(hdr, "##unused=unformatted text 2");
    bcf_hdr_append(hdr, "##contig=<ID=Unused,length=62435964>");
    bcf_hdr_append(hdr, "##source=myImputationProgramV3.1");
    bcf_hdr_append(hdr, "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta");
    bcf_hdr_append(hdr, "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>");
    bcf_hdr_append(hdr, "##phasing=partial");
    bcf_hdr_append(hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
    bcf_hdr_append(hdr, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
    bcf_hdr_append(hdr, "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
    bcf_hdr_append(hdr, "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=q10,Description=\"Quality below 10\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=TS,Number=1,Type=String,Description=\"Test String\">");

    bcf_hdr_add_sample(hdr, "NA00001");
    bcf_hdr_add_sample(hdr, "NA00002");
    bcf_hdr_add_sample(hdr, "NA00003");
    bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
    bcf_hdr_write(fp, hdr);


    // Add a record
    // 20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
    // .. CHROM
    rec->rid = bcf_hdr_name2id(hdr, "20");
    // .. POS
    rec->pos = 14369;
    // .. ID
    bcf_update_id(hdr, rec, "rs6054257");
    // .. REF and ALT
    bcf_update_alleles_str(hdr, rec, "G,A");
    // .. QUAL
    rec->qual = 29;
    // .. FILTER
    int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
    bcf_update_filter(hdr, rec, &tmpi, 1);
    // .. INFO
    tmpi = 3;
    bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
    tmpi = 14;
    bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
    float tmpf = 0.5;
    bcf_update_info_float(hdr, rec, "AF", &tmpf, 1);
    bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
    bcf_update_info_flag(hdr, rec, "H2", NULL, 1);
    // .. FORMAT
    int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
    tmpia[0] = bcf_gt_phased(0);
    tmpia[1] = bcf_gt_phased(0);
    tmpia[2] = bcf_gt_phased(1);
    tmpia[3] = bcf_gt_phased(0);
    tmpia[4] = bcf_gt_unphased(1);
    tmpia[5] = bcf_gt_unphased(1);
    bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
    tmpia[0] = 48;
    tmpia[1] = 48;
    tmpia[2] = 43;
    bcf_update_format_int32(hdr, rec, "GQ", tmpia, bcf_hdr_nsamples(hdr));
    tmpia[0] = 1;
    tmpia[1] = 8;
    tmpia[2] = 5;
    bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr));
    tmpia[0] = 51;
    tmpia[1] = 51;
    tmpia[2] = 51;
    tmpia[3] = 51;
    tmpia[4] = bcf_int32_missing;
    tmpia[5] = bcf_int32_missing;
    bcf_update_format_int32(hdr, rec, "HQ", tmpia, bcf_hdr_nsamples(hdr)*2);
    char *tmp_str[] = {"String1","SomeOtherString2","YetAnotherString3"};
    bcf_update_format_string(hdr, rec, "TS", (const char**)tmp_str, 3);
    bcf_write1(fp, hdr, rec);

    // 20     1110696 . A      G,T     67   .   NS=2;DP=10;AF=0.333,.;AA=T;DB GT 2 1   ./.
    bcf_clear1(rec);
    rec->rid = bcf_hdr_name2id(hdr, "20");
    rec->pos = 1110695;
    bcf_update_alleles_str(hdr, rec, "A,G,T");
    rec->qual = 67;
    tmpi = 2;
    bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
    tmpi = 10;
    bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
    float *tmpfa = (float*)malloc(2*sizeof(float));
    tmpfa[0] = 0.333;
    bcf_float_set_missing(tmpfa[1]);
    bcf_update_info_float(hdr, rec, "AF", tmpfa, 2);
    bcf_update_info_string(hdr, rec, "AA", "T");
    bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
    tmpia[0] = bcf_gt_phased(2);
    tmpia[1] = bcf_int32_vector_end;
    tmpia[2] = bcf_gt_phased(1);
    tmpia[3] = bcf_int32_vector_end;
    tmpia[4] = bcf_gt_missing;
    tmpia[5] = bcf_gt_missing;
    bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
    bcf_write1(fp, hdr, rec);

    free(tmpia);
    free(tmpfa);

    // Clean
    free(str.s);
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
}

void bcf_to_vcf(char *fname)
{
    htsFile *fp    = hts_open(fname,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec    = bcf_init1();

    char *gz_fname = (char*) malloc(strlen(fname)+4);
    snprintf(gz_fname,strlen(fname)+4,"%s.gz",fname);
    htsFile *out   = hts_open(gz_fname,"wg");

    bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
    bcf_hdr_remove(hdr_out,BCF_HL_STR,"unused");
    bcf_hdr_remove(hdr_out,BCF_HL_GEN,"unused");
    bcf_hdr_remove(hdr_out,BCF_HL_FLT,"Flt");
    bcf_hdr_remove(hdr_out,BCF_HL_INFO,"UI");
    bcf_hdr_remove(hdr_out,BCF_HL_FMT,"UF");
    bcf_hdr_remove(hdr_out,BCF_HL_CTG,"Unused");
    bcf_hdr_write(out, hdr_out);

    while ( bcf_read1(fp, hdr, rec)>=0 )
    {
        bcf_write1(out, hdr_out, rec);

        // Test problems caused by bcf1_sync: the data block
        // may be realloced, also the unpacked structures must
        // get updated.
        bcf_unpack(rec, BCF_UN_STR);
        bcf_update_id(hdr, rec, 0);
        bcf_update_format_int32(hdr, rec, "GQ", NULL, 0);

        bcf1_t *dup = bcf_dup(rec);     // force bcf1_sync call
        bcf_write1(out, hdr_out, dup);
        bcf_destroy1(dup);

        bcf_update_alleles_str(hdr_out, rec, "G,A");
        int32_t tmpi = 99;
        bcf_update_info_int32(hdr_out, rec, "DP", &tmpi, 1);
        int32_t tmpia[] = {9,9,9};
        bcf_update_format_int32(hdr_out, rec, "DP", tmpia, 3);

        bcf_write1(out, hdr_out, rec);
    }

    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    bcf_hdr_destroy(hdr_out);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
    if ( (ret=hts_close(out)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",gz_fname,ret);
        exit(ret);
    }


    // read gzip, write stdout
    htsFile *gz_in = hts_open(gz_fname, "r");
    if ( !gz_in )
    {
        fprintf(stderr,"Could not read: %s\n", gz_fname);
        exit(1);
    }

    kstring_t line = {0,0,0};
    while ( hts_getline(gz_in, KS_SEP_LINE, &line)>0 )
    {
        kputc('\n',&line);
        fwrite(line.s,1,line.l,stdout);
    }

    if ( (ret=hts_close(gz_in)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",gz_fname,ret);
        exit(ret);
    }
    free(line.s);
    free(gz_fname);
}

void iterator(const char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    hts_idx_t *idx;
    hts_itr_t *iter;

    bcf_index_build(fname, 0);
    idx = bcf_index_load(fname);

    iter = bcf_itr_queryi(idx, bcf_hdr_name2id(hdr, "20"), 1110600, 1110800);
    bcf_itr_destroy(iter);

    iter = bcf_itr_querys(idx, hdr, "20:1110600-1110800");
    bcf_itr_destroy(iter);

    hts_idx_destroy(idx);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
}

void test_get_info_values(const char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *line = bcf_init();

    while (bcf_read(fp, hdr, line) == 0)
    {
        float *afs = 0;
        int count = 0;
        int ret = bcf_get_info_float(hdr, line, "AF", &afs, &count);

        if (line->pos == 14369)
        {
            if (ret != 1 || afs[0] != 0.5f)
            {
                fprintf(stderr, "AF on position 14370 should be 0.5\n");
                exit(-1);
            }
        }
        else
        {
            if (ret != 2 || afs[0] != 0.333f || !bcf_float_is_missing(afs[1]))
            {
                fprintf(stderr, "AF on position 1110696 should be 0.333, missing\n");
                exit(-1);
            }
        }

        free(afs);
    }

    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
}

void write_format_values(const char *fname)
{
    // Init
    htsFile *fp = hts_open(fname, "wb");
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    bcf1_t *rec = bcf_init1();

    // Create VCF header
    bcf_hdr_append(hdr, "##contig=<ID=1>");
    bcf_hdr_append(hdr, "##FORMAT=<ID=TF,Number=1,Type=Float,Description=\"Test Float\">");
    bcf_hdr_add_sample(hdr, "S");
    bcf_hdr_add_sample(hdr, NULL); // to update internal structures
    bcf_hdr_write(fp, hdr);

    // Add a record
    // .. FORMAT
    float test[4];
    bcf_float_set_missing(test[0]);
    test[1] = 47.11f;
    bcf_float_set_vector_end(test[2]);
    test[3] = -1.2e-13;
    bcf_update_format_float(hdr, rec, "TF", test, 4);
    bcf_write1(fp, hdr, rec);

    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ((ret = hts_close(fp)))
    {
        fprintf(stderr, "hts_close(%s): non-zero status %d\n", fname, ret);
        exit(ret);
    }
}

void check_format_values(const char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *line = bcf_init();

    while (bcf_read(fp, hdr, line) == 0)
    {
        float *values = 0;
        int count = 0;
        int ret = bcf_get_format_float(hdr, line, "TF", &values, &count);

        // NOTE the return value from bcf_get_format_float is different from
        // bcf_get_info_float in the sense that vector-end markers also count.
        if (ret != 4 ||
            count < ret ||
            !bcf_float_is_missing(values[0]) ||
            values[1] != 47.11f ||
            !bcf_float_is_vector_end(values[2]) ||
            !bcf_float_is_vector_end(values[3]))
        {
            fprintf(stderr, "bcf_get_format_float didn't produce the expected output.\n");
            exit(-1);
        }

        free(values);
    }

    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
}

void test_get_format_values(const char *fname)
{
    write_format_values(fname);
    check_format_values(fname);
}

int main(int argc, char **argv)
{
    char *fname = argc>1 ? argv[1] : "rmme.bcf";

    // format test. quiet unless there's a failure
    test_get_format_values(fname);

    // main test. writes to stdout
    write_bcf(fname);
    bcf_to_vcf(fname);
    iterator(fname);

    // additional tests. quiet unless there's a failure.
    test_get_info_values(fname);

    return 0;
}
