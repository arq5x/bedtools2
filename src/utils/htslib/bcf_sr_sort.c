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

#include <config.h>

#include <strings.h>

#include "bcf_sr_sort.h"
#include "htslib/khash_str2int.h"

#define SR_REF   1
#define SR_SNP   2
#define SR_INDEL 4
#define SR_OTHER 8
#define SR_SCORE(srt,a,b) (srt)->score[((a)<<4)|(b)]

// Resize a bit set.
static inline kbitset_t *kbs_resize(kbitset_t *bs, size_t ni)
{
    if ( !bs ) return kbs_init(ni);
    size_t n = (ni + KBS_ELTBITS-1) / KBS_ELTBITS;
    if ( n==bs->n ) return bs;

    bs = (kbitset_t *) realloc(bs, sizeof(kbitset_t) + n * sizeof(unsigned long));
    if ( bs==NULL ) return NULL;
    if ( n > bs->n )
        memset(bs->b + bs->n, 0, (n - bs->n) * sizeof (unsigned long));
    bs->n = n;
    bs->b[n] = ~0UL;
    return bs;
}

// Logical AND
static inline int kbs_logical_and(kbitset_t *bs1, kbitset_t *bs2)
{
    // General case, bitsets of unequal size:
    //  int i, n = bs1->n < bs2->n ? bs1->n : bs2->n;
    int i, n = bs1->n;

    for (i=0; i<n; i++) if ( bs1->b[i] & bs2->b[i] ) return 1;
    return 0;
}

// Bitwise OR, dst will be modified, src will be left unchanged
static inline void kbs_bitwise_or(kbitset_t *dst, kbitset_t *src)
{
    int i;
    for (i=0; i<dst->n; i++) dst->b[i] |= src->b[i];
}


static void bcf_sr_init_scores(sr_sort_t *srt)
{
    int i,jbit,kbit;

    // lower number = lower priority, zero means forbidden

    if ( srt->pair & BCF_SR_PAIR_ANY ) srt->pair |= (BCF_SR_PAIR_SNPS | BCF_SR_PAIR_INDELS | BCF_SR_PAIR_SNP_REF | BCF_SR_PAIR_INDEL_REF);
    if ( srt->pair & BCF_SR_PAIR_SNPS ) SR_SCORE(srt,SR_SNP,SR_SNP) = 3;
    if ( srt->pair & BCF_SR_PAIR_INDELS ) SR_SCORE(srt,SR_INDEL,SR_INDEL) = 3;
    if ( srt->pair & BCF_SR_PAIR_SNP_REF )
    {
        SR_SCORE(srt,SR_SNP,SR_REF) = 2;
        SR_SCORE(srt,SR_REF,SR_SNP) = 2;
    }
    if ( srt->pair & BCF_SR_PAIR_INDEL_REF )
    {
        SR_SCORE(srt,SR_INDEL,SR_REF) = 2;
        SR_SCORE(srt,SR_REF,SR_INDEL) = 2;
    }
    if ( srt->pair & BCF_SR_PAIR_ANY )
    {
        for (i=0; i<256; i++)
            if ( !srt->score[i] ) srt->score[i] = 1;
    }

    // set all combinations
    for (i=0; i<256; i++)
    {
        if ( srt->score[i] ) continue;      // already set
        int max = 0;
        for (jbit=0; jbit<4; jbit++)        // high bits
        {
            int j = 1<<jbit;
            if ( !(i & (j<<4)) ) continue;
            for (kbit=0; kbit<4; kbit++)    // low bits
            {
                int k = 1<<kbit;
                if ( !(i & k) ) continue;
                if ( max < SR_SCORE(srt,j,k) ) max = SR_SCORE(srt,j,k);
            }
        }
        srt->score[i] = max;
    }
}
static int multi_is_exact(var_t *avar, var_t *bvar)
{
    if ( avar->nalt != bvar->nalt ) return 0;

    int alen = strlen(avar->str);
    int blen = strlen(bvar->str);
    if ( alen != blen ) return 0;

    char *abeg = avar->str;
    while ( *abeg )
    {
        char *aend = abeg;
        while ( *aend && *aend!=',' ) aend++;

        char *bbeg = bvar->str;
        while ( *bbeg )
        {
            char *bend = bbeg;
            while ( *bend && *bend!=',' ) bend++;
            if ( bend - bbeg == aend - abeg && !strncasecmp(abeg,bbeg,bend-bbeg) ) break;
            bbeg = *bend ? bend+1 : bend;
        }
        if ( !*bbeg ) return 0;

        abeg = *aend ? aend+1 : aend;
    }
    return 1;
}
static int multi_is_subset(var_t *avar, var_t *bvar)
{
    char *abeg = avar->str;
    while ( *abeg )
    {
        char *aend = abeg;
        while ( *aend && *aend!=',' ) aend++;

        char *bbeg = bvar->str;
        while ( *bbeg )
        {
            char *bend = bbeg;
            while ( *bend && *bend!=',' ) bend++;
            if ( bend - bbeg == aend - abeg && !strncasecmp(abeg,bbeg,bend-bbeg) ) return 1;
            bbeg = *bend ? bend+1 : bend;
        }
        abeg = *aend ? aend+1 : aend;
    }
    return 0;
}
int32_t pairing_score(sr_sort_t *srt, int ivset, int jvset)
{
    varset_t *iv = &srt->vset[ivset];
    varset_t *jv = &srt->vset[jvset];

    // Restrictive logic: the strictest type from a group is selected,
    // so that, for example, snp+ref does not lead to the inclusion of an indel
    int i,j;
    uint32_t min = UINT32_MAX;
    for (i=0; i<iv->nvar; i++)
    {
        var_t *ivar = &srt->var[iv->var[i]];
        for (j=0; j<jv->nvar; j++)
        {
            var_t *jvar = &srt->var[jv->var[j]];
            if ( srt->pair & BCF_SR_PAIR_EXACT )
            {
                if ( ivar->type != jvar->type ) continue;
                if ( !strcmp(ivar->str,jvar->str) ) return UINT32_MAX;  // exact match, best possibility
                if ( multi_is_exact(ivar,jvar) ) return UINT32_MAX; // identical alleles
                continue;
            }
            if ( ivar->type==jvar->type && !strcmp(ivar->str,jvar->str) ) return UINT32_MAX;  // exact match, best possibility
            if ( ivar->type & jvar->type && multi_is_subset(ivar,jvar) ) return UINT32_MAX; // one of the alleles is identical

            uint32_t score = SR_SCORE(srt,ivar->type,jvar->type);
            if ( !score ) return 0;     // some of the varsets in the two groups are not compatible, will not pair
            if ( min>score ) min = score;
        }
    }
    if ( srt->pair & BCF_SR_PAIR_EXACT ) return 0;

    assert( min!=UINT32_MAX );

    uint32_t cnt = 0;
    for (i=0; i<iv->nvar; i++) cnt += srt->var[iv->var[i]].nvcf;
    for (j=0; j<jv->nvar; j++) cnt += srt->var[jv->var[j]].nvcf;

    return (1<<(28+min)) + cnt;
}
void remove_vset(sr_sort_t *srt, int jvset)
{
    if ( jvset+1 < srt->nvset )
    {
        varset_t tmp = srt->vset[jvset];
        memmove(&srt->vset[jvset], &srt->vset[jvset+1], sizeof(varset_t)*(srt->nvset - jvset - 1));
        srt->vset[srt->nvset-1] = tmp;

        int *jmat = srt->pmat + jvset*srt->ngrp;
        memmove(jmat, &jmat[srt->ngrp],sizeof(int)*(srt->nvset - jvset - 1)*srt->ngrp);

        memmove(&srt->cnt[jvset], &srt->cnt[jvset+1], sizeof(int)*(srt->nvset - jvset - 1));
    }
    srt->nvset--;
}
int merge_vsets(sr_sort_t *srt, int ivset, int jvset)
{
    int i,j;
    if ( ivset > jvset ) { i = ivset; ivset = jvset; jvset = i; }

    varset_t *iv = &srt->vset[ivset];
    varset_t *jv = &srt->vset[jvset];

    kbs_bitwise_or(iv->mask,jv->mask);

    i = iv->nvar;
    iv->nvar += jv->nvar;
    hts_expand(int, iv->nvar, iv->mvar, iv->var);
    for (j=0; j<jv->nvar; j++,i++) iv->var[i] = jv->var[j];

    int *imat = srt->pmat + ivset*srt->ngrp;
    int *jmat = srt->pmat + jvset*srt->ngrp;
    for (i=0; i<srt->ngrp; i++) imat[i] += jmat[i];
    srt->cnt[ivset] += srt->cnt[jvset];

    remove_vset(srt, jvset);

    return ivset;
}
void push_vset(sr_sort_t *srt, int ivset)
{
    varset_t *iv = &srt->vset[ivset];
    int i,j;
    for (i=0; i<srt->sr->nreaders; i++)
    {
        vcf_buf_t *buf = &srt->vcf_buf[i];
        buf->nrec++;
        hts_expand(bcf1_t*,buf->nrec,buf->mrec,buf->rec);
        buf->rec[buf->nrec-1] = NULL;
    }
    for (i=0; i<iv->nvar; i++)
    {
        var_t *var = &srt->var[ iv->var[i] ];
        for (j=0; j<var->nvcf; j++)
        {
            int jvcf = var->vcf[j];
            vcf_buf_t *buf = &srt->vcf_buf[jvcf];
            buf->rec[buf->nrec-1] = var->rec[j];
        }
    }
    remove_vset(srt, ivset);
}

static int cmpstringp(const void *p1, const void *p2)
{
    return strcmp(* (char * const *) p1, * (char * const *) p2);
}

#if DEBUG_VSETS
void debug_vsets(sr_sort_t *srt)
{
    int i,j,k;
    for (i=0; i<srt->nvset; i++)
    {
        fprintf(stderr,"dbg_vset %d:", i);
        for (j=0; j<srt->vset[i].mask->n; j++) fprintf(stderr,"%c%lu",j==0?' ':':',srt->vset[i].mask->b[j]);
        fprintf(stderr,"\t");
        for (j=0; j<srt->vset[i].nvar; j++)
        {
            var_t *var = &srt->var[srt->vset[i].var[j]];
            fprintf(stderr,"\t%s",var->str);
            for (k=0; k<var->nvcf; k++)
                fprintf(stderr,"%c%d", k==0?':':',',var->vcf[k]);
        }
        fprintf(stderr,"\n");
    }
}
#endif

#if DEBUG_VBUF
void debug_vbuf(sr_sort_t *srt)
{
    int i, j;
    for (j=0; j<srt->vcf_buf[0].nrec; j++)
    {
        fprintf(stderr,"dbg_vbuf %d:\t", j);
        for (i=0; i<srt->sr->nreaders; i++)
        {
            vcf_buf_t *buf = &srt->vcf_buf[i];
            fprintf(stderr,"\t%d", buf->rec[j] ? buf->rec[j]->pos+1 : 0);
        }
        fprintf(stderr,"\n");
    }
}
#endif

char *grp_create_key(sr_sort_t *srt)
{
    if ( !srt->str.l ) return strdup("");
    int i;
    hts_expand(char*,srt->noff,srt->mcharp,srt->charp);
    for (i=0; i<srt->noff; i++)
    {
        srt->charp[i] = srt->str.s + srt->off[i];
        if ( i>0 ) srt->charp[i][-1] = 0;
    }
    qsort(srt->charp, srt->noff, sizeof(*srt->charp), cmpstringp);
    char *ret = (char*) malloc(srt->str.l + 1), *ptr = ret;
    for (i=0; i<srt->noff; i++)
    {
        int len = strlen(srt->charp[i]);
        memcpy(ptr, srt->charp[i], len);
        ptr += len + 1;
        ptr[-1] = i+1==srt->noff ? 0 : ';';
    }
    return ret;
}
int bcf_sr_sort_set_active(sr_sort_t *srt, int idx)
{
    hts_expand(int,idx+1,srt->mactive,srt->active);
    srt->nactive = 1;
    srt->active[srt->nactive - 1] = idx;
    return 0;
}
int bcf_sr_sort_add_active(sr_sort_t *srt, int idx)
{
    hts_expand(int,idx+1,srt->mactive,srt->active);
    srt->nactive++;
    srt->active[srt->nactive - 1] = idx;
    return 0;
}
static void bcf_sr_sort_set(bcf_srs_t *readers, sr_sort_t *srt, const char *chr, int min_pos)
{
    if ( !srt->grp_str2int )
    {
        // first time here, initialize
        if ( !srt->pair )
        {
            if ( readers->collapse==COLLAPSE_NONE ) readers->collapse = BCF_SR_PAIR_EXACT;
            bcf_sr_set_opt(readers, BCF_SR_PAIR_LOGIC, readers->collapse);
        }
        bcf_sr_init_scores(srt);
        srt->grp_str2int = khash_str2int_init();
        srt->var_str2int = khash_str2int_init();
    }
    int k;
    khash_t(str2int) *hash;
    hash = srt->grp_str2int;
    for (k=0; k < kh_end(hash); k++)
        if ( kh_exist(hash,k) ) free((char*)kh_key(hash,k));
    hash = srt->var_str2int;
    for (k=0; k < kh_end(hash); k++)
        if ( kh_exist(hash,k) ) free((char*)kh_key(hash,k));
    kh_clear(str2int, srt->grp_str2int);
    kh_clear(str2int, srt->var_str2int);
    srt->ngrp = srt->nvar = srt->nvset = 0;

    grp_t grp;
    memset(&grp,0,sizeof(grp_t));

    // group VCFs into groups, each with a unique combination of variants in the duplicate lines
    int ireader,ivar,irec,igrp,ivset,iact;
    for (ireader=0; ireader<readers->nreaders; ireader++) srt->vcf_buf[ireader].nrec = 0;
    for (iact=0; iact<srt->nactive; iact++)
    {
        ireader = srt->active[iact];
        bcf_sr_t *reader = &readers->readers[ireader];
        int rid   = bcf_hdr_name2id(reader->header, chr);
        grp.nvar  = 0;
        hts_expand(int,reader->nbuffer,srt->moff,srt->off);
        srt->noff  = 0;
        srt->str.l = 0;
        for (irec=1; irec<=reader->nbuffer; irec++)
        {
            bcf1_t *line = reader->buffer[irec];
            if ( line->rid!=rid || line->pos!=min_pos ) break;

            if ( srt->str.l ) kputc(';',&srt->str);
            srt->off[srt->noff++] = srt->str.l;
            size_t beg = srt->str.l;
            for (ivar=1; ivar<line->n_allele; ivar++)
            {
                if ( ivar>1 ) kputc(',',&srt->str);
                kputs(line->d.allele[0],&srt->str);
                kputc('>',&srt->str);
                kputs(line->d.allele[ivar],&srt->str);
            }
            if ( line->n_allele==1 )
            {
                kputs(line->d.allele[0],&srt->str);
                kputsn(">.",2,&srt->str);
            }

            // Create new variant or attach to existing one. But careful, there can be duplicate
            // records with the same POS,REF,ALT (e.g. in dbSNP-b142)
            char *var_str = beg + srt->str.s;
            int ret, var_idx = 0, var_end = srt->str.l;
            while ( 1 )
            {
                ret = khash_str2int_get(srt->var_str2int, var_str, &ivar);
                if ( ret==-1 ) break;

                var_t *var = &srt->var[ivar];
                if ( var->vcf[var->nvcf-1] != ireader ) break;

                srt->str.l = var_end;
                kputw(var_idx, &srt->str);
                var_str = beg + srt->str.s;
                var_idx++;
            }
            if ( ret==-1 )
            {
                ivar = srt->nvar++;
                hts_expand0(var_t,srt->nvar,srt->mvar,srt->var);
                srt->var[ivar].nvcf = 0;
                khash_str2int_set(srt->var_str2int, strdup(var_str), ivar);
                free(srt->var[ivar].str);   // possible left-over from the previous position
            }
            var_t *var = &srt->var[ivar];
            var->nalt = line->n_allele - 1;
            var->type = bcf_get_variant_types(line);
            srt->str.s[var_end] = 0;
            if ( ret==-1 )
                var->str = strdup(var_str);

            int mvcf = var->mvcf;
            var->nvcf++;
            hts_expand0(int*, var->nvcf, var->mvcf, var->vcf);
            if ( mvcf != var->mvcf ) var->rec = (bcf1_t **) realloc(var->rec,sizeof(bcf1_t*)*var->mvcf);
            var->vcf[var->nvcf-1] = ireader;
            var->rec[var->nvcf-1] = line;

            grp.nvar++;
            hts_expand(var_t,grp.nvar,grp.mvar,grp.var);
            grp.var[grp.nvar-1] = ivar;
        }
        char *grp_key = grp_create_key(srt);
        int ret = khash_str2int_get(srt->grp_str2int, grp_key, &igrp);
        if ( ret==-1 )
        {
            igrp = srt->ngrp++;
            hts_expand0(grp_t, srt->ngrp, srt->mgrp, srt->grp);
            free(srt->grp[igrp].var);
            srt->grp[igrp] = grp;
            srt->grp[igrp].key = grp_key;
            khash_str2int_set(srt->grp_str2int, grp_key, igrp);
            memset(&grp,0,sizeof(grp_t));
        }
        else
            free(grp_key);
        srt->grp[igrp].nvcf++;
    }
    free(grp.var);

    // initialize bitmask - which groups is the variant present in
    for (ivar=0; ivar<srt->nvar; ivar++)
    {
        srt->var[ivar].mask = kbs_resize(srt->var[ivar].mask, srt->ngrp);
        kbs_clear(srt->var[ivar].mask);
    }
    for (igrp=0; igrp<srt->ngrp; igrp++)
    {
        for (ivar=0; ivar<srt->grp[igrp].nvar; ivar++)
        {
            int i = srt->grp[igrp].var[ivar];
            kbs_insert(srt->var[i].mask, igrp);
        }
    }

    // create the initial list of variant sets
    for (ivar=0; ivar<srt->nvar; ivar++)
    {
        ivset = srt->nvset++;
        hts_expand0(varset_t, srt->nvset, srt->mvset, srt->vset);

        varset_t *vset = &srt->vset[ivset];
        vset->nvar = 1;
        hts_expand0(var_t, vset->nvar, vset->mvar, vset->var);
        vset->var[vset->nvar-1] = ivar;
        var_t *var  = &srt->var[ivar];
        vset->cnt   = var->nvcf;
        vset->mask  = kbs_resize(vset->mask, srt->ngrp);
        kbs_clear(vset->mask);
        kbs_bitwise_or(vset->mask, var->mask);

        int type = 0;
        if ( var->type==VCF_REF ) type |= SR_REF;
        else
        {
            if ( var->type & VCF_SNP ) type |= SR_SNP;
            if ( var->type & VCF_MNP ) type |= SR_SNP;
            if ( var->type & VCF_INDEL ) type |= SR_INDEL;
            if ( var->type & VCF_OTHER ) type |= SR_OTHER;
        }
        var->type = type;
    }
#if DEBUG_VSETS
    debug_vsets(srt);
#endif

    // initialize the pairing matrix
    hts_expand(int, srt->ngrp*srt->nvset, srt->mpmat, srt->pmat);
    hts_expand(int, srt->nvset, srt->mcnt, srt->cnt);
    memset(srt->pmat, 0, sizeof(*srt->pmat)*srt->ngrp*srt->nvset);
    for (ivset=0; ivset<srt->nvset; ivset++)
    {
        varset_t *vset = &srt->vset[ivset];
        for (igrp=0; igrp<srt->ngrp; igrp++) srt->pmat[ivset*srt->ngrp+igrp] = 0;
        srt->cnt[ivset] = vset->cnt;
    }

    // pair the lines
    while ( srt->nvset )
    {
#if DEBUG_VSETS
    fprintf(stderr,"\n");
    debug_vsets(srt);
#endif

        int imax = 0;
        for (ivset=1; ivset<srt->nvset; ivset++)
            if ( srt->cnt[imax] < srt->cnt[ivset] ) imax = ivset;

        int ipair = -1;
        uint32_t max_score = 0;
        for (ivset=0; ivset<srt->nvset; ivset++)
        {
            if ( kbs_logical_and(srt->vset[imax].mask,srt->vset[ivset].mask) ) continue;   // cannot be merged
            uint32_t score = pairing_score(srt, imax, ivset);
            // fprintf(stderr,"score: %d %d, logic=%d \t..\t %u\n", imax,ivset,srt->pair,score);
            if ( max_score < score ) { max_score = score; ipair = ivset; }
        }

        // merge rows creating a new variant set this way
        if ( ipair!=-1 && ipair!=imax )
        {
            imax = merge_vsets(srt, imax, ipair);
            continue;
        }

        push_vset(srt, imax);
    }

    srt->chr = chr;
    srt->pos = min_pos;
}

int bcf_sr_sort_next(bcf_srs_t *readers, sr_sort_t *srt, const char *chr, int min_pos)
{
    int i,j;
    assert( srt->nactive>0 );

    if ( srt->nsr != readers->nreaders )
    {
        srt->sr = readers;
        if ( srt->nsr < readers->nreaders )
        {
            srt->vcf_buf = (vcf_buf_t*) realloc(srt->vcf_buf,readers->nreaders*sizeof(vcf_buf_t));
            memset(srt->vcf_buf + srt->nsr, 0, sizeof(vcf_buf_t)*(readers->nreaders - srt->nsr));
            if ( srt->msr < srt->nsr ) srt->msr = srt->nsr;
        }
        srt->nsr = readers->nreaders;
        srt->chr = NULL;
    }
    if ( srt->nactive == 1 )
    {
        if ( readers->nreaders>1 )
            memset(readers->has_line, 0, readers->nreaders*sizeof(*readers->has_line));
        bcf_sr_t *reader = &readers->readers[srt->active[0]];
        assert( reader->buffer[1]->pos==min_pos );
        bcf1_t *tmp = reader->buffer[0];
        for (j=1; j<=reader->nbuffer; j++) reader->buffer[j-1] = reader->buffer[j];
        reader->buffer[ reader->nbuffer ] = tmp;
        reader->nbuffer--;
        readers->has_line[srt->active[0]] = 1;
        return 1;
    }
    if ( !srt->chr || srt->pos!=min_pos || strcmp(srt->chr,chr) ) bcf_sr_sort_set(readers, srt, chr, min_pos);

    if ( !srt->vcf_buf[0].nrec ) return 0;

#if DEBUG_VBUF
    debug_vbuf(srt);
#endif

    int nret = 0;
    for (i=0; i<srt->sr->nreaders; i++)
    {
        vcf_buf_t *buf = &srt->vcf_buf[i];

        if ( buf->rec[0] )
        {
            bcf_sr_t *reader = &srt->sr->readers[i];
            for (j=1; j<=reader->nbuffer; j++)
                if ( reader->buffer[j] == buf->rec[0] ) break;

            assert( j<=reader->nbuffer );

            bcf1_t *tmp = reader->buffer[0];
            reader->buffer[0] = reader->buffer[j++];
            for (; j<=reader->nbuffer; j++) reader->buffer[j-1] = reader->buffer[j];
            reader->buffer[ reader->nbuffer ] = tmp;
            reader->nbuffer--;

            nret++;
            srt->sr->has_line[i] = 1;
        }
        else
            srt->sr->has_line[i] = 0;

        buf->nrec--;
        if ( buf->nrec > 0 )
            memmove(buf->rec, &buf->rec[1], buf->nrec*sizeof(bcf1_t*));
    }
    return nret;
}
void bcf_sr_sort_remove_reader(bcf_srs_t *readers, sr_sort_t *srt, int i)
{
    //vcf_buf is allocated only in bcf_sr_sort_next
    //So, a call to bcf_sr_add_reader() followed immediately by bcf_sr_remove_reader()
    //would cause the program to crash in this segment
    if (srt->vcf_buf)
    {
        free(srt->vcf_buf[i].rec);
        if ( i+1 < srt->nsr )
            memmove(&srt->vcf_buf[i], &srt->vcf_buf[i+1], (srt->nsr - i - 1)*sizeof(vcf_buf_t));
        memset(srt->vcf_buf + srt->nsr - 1, 0, sizeof(vcf_buf_t));
    }
}
sr_sort_t *bcf_sr_sort_init(sr_sort_t *srt)
{
    if ( !srt ) return calloc(1,sizeof(sr_sort_t));
    memset(srt,0,sizeof(sr_sort_t));
    return srt;
}
void bcf_sr_sort_reset(sr_sort_t *srt)
{
    srt->chr = NULL;
}
void bcf_sr_sort_destroy(sr_sort_t *srt)
{
    free(srt->active);
    if ( srt->var_str2int ) khash_str2int_destroy_free(srt->var_str2int);
    if ( srt->grp_str2int ) khash_str2int_destroy_free(srt->grp_str2int);
    int i;
    for (i=0; i<srt->nsr; i++) free(srt->vcf_buf[i].rec);
    free(srt->vcf_buf);
    for (i=0; i<srt->mvar; i++)
    {
        free(srt->var[i].str);
        free(srt->var[i].vcf);
        free(srt->var[i].rec);
        kbs_destroy(srt->var[i].mask);
    }
    free(srt->var);
    for (i=0; i<srt->mgrp; i++)
        free(srt->grp[i].var);
    free(srt->grp);
    for (i=0; i<srt->mvset; i++)
    {
        kbs_destroy(srt->vset[i].mask);
        free(srt->vset[i].var);
    }
    free(srt->vset);
    free(srt->str.s);
    free(srt->off);
    free(srt->charp);
    free(srt->cnt);
    free(srt->pmat);
    memset(srt,0,sizeof(*srt));
}

