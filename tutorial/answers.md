% bedtools Tutorial
% Aaron Quinlan


Puzzles to help teach you more bedtools.
========================================

1. Create a BED file representing all of the intervals in the genome
that are NOT exonic and are not Promoters (based on the promoters in the hESC file).

Answer:
     
    grep Promoter hesc.chromHmm.bed > hesc.promoters.bed
    
    cat exons.bed hesc.promoters.bed | sort -k1,1 -k2,2n | exons.and.promoters.bed

    bedtools complement -i exons.and.promoters.bed -g genome.txt > notexonsorpromoters.bed


2. What is the average distance from GWAS SNPs to the closest exon? (Hint - have a look at the [closest](http://bedtools.readthedocs.org/en/latest/content/tools/closest.html) tool.)

Answer:

    bedtools closest -a gwas.bed -b exons.bed -d | head
    chr1	1005805	1005806	rs3934834	chr1	1007125	1007955	NM_001205252_exon_0_0_chr1_1007126_r	0	-	1320
    chr1	1079197	1079198	rs11260603	chr1	1078118	1079434	NR_038869_exon_2_0_chr1_1078119_f	0	+	0
    chr1	1247493	1247494	rs12103	chr1	1247397	1247527	NM_001256456_exon_1_0_chr1_1247398_r	0	-	0
    chr1	1247493	1247494	rs12103	chr1	1247397	1247527	NM_001256460_exon_1_0_chr1_1247398_r	0	-	0
    chr1	1247493	1247494	rs12103	chr1	1247397	1247527	NM_001256462_exon_1_0_chr1_1247398_r	0	-	0
    chr1	1247493	1247494	rs12103	chr1	1247397	1247527	NM_001256463_exon_1_0_chr1_1247398_r	0	-	0
    chr1	1247493	1247494	rs12103	chr1	1247397	1247527	NM_017871_exon_1_0_chr1_1247398_r	0	-	0
    chr1	2069171	2069172	rs425277	chr1	2066700	2066786	NM_001033581_exon_1_0_chr1_2066701_f	0	+	2386
    chr1	2069171	2069172	rs425277	chr1	2066700	2066786	NM_001033582_exon_1_0_chr1_2066701_f	0	+	2386
    chr1	2069171	2069172	rs425277	chr1	2066700	2066786	NM_001242874_exon_1_0_chr1_2066701_f	0	+	2386

    bedtools closest -a gwas.bed -b exons.bed -d \
      | awk '{ sum += $11 } END { if (NR > 0) print sum / NR }'
    46713.1

3. Count how many exons occur in each 500kb interval ("window") in the human genome. (Hint - have a look at the `makewindows` tool.)

Answer:

    bedtools makewindows -g genome.txt -w 500000 > genome.windows.bed
    bedtools intersect -a genome.windows.bed -b exons.bed -c > genome.windows.exoncount.bedg

or...

	bedtools makewindows -g genome.txt -w 500000 \
      | bedtools intersect -a - -b exons.bed -c \
    > genome.windows.exoncount.bedg

4. Are there any exons that are completely overlapped by an enhancer? If so, how many?

Answer:

    bedtools intersect -a exons.bed \
                       -b <(grep Enhancer hesc.chromHmm.bed) \
                       -wa -wb -f 1.0 \
    | head
    chr1	948846	948956	NM_005101_exon_0_0_chr1_948847_f	0	+	chr1	948337	949337	4_Strong_Enhancer
    chr1	1051439	1051736	NM_017891_exon_9_0_chr1_1051440_r	0	-	chr1	1051337	1051737	6_Weak_Enhancer
    chr1	1109285	1109306	NM_001130045_exon_0_0_chr1_1109286_f	0	+	chr1	1108537	1109537	6_Weak_Enhancer
    chr1	1109803	1109869	NM_001130045_exon_2_0_chr1_1109804_f	0	+	chr1	1109737	1109937	6_Weak_Enhancer
    chr1	1219357	1219470	NM_001130413_exon_4_0_chr1_1219358_f	0	+	chr1	1219137	1220137	7_Weak_Enhancer
    chr1	1219357	1219470	NR_037668_exon_4_0_chr1_1219358_f	0	+	chr1	1219137	1220137	7_Weak_Enhancer
    chr1	1229202	1229313	NM_030649_exon_1_0_chr1_1229203_r	0	-	chr1	1228937	1229937	6_Weak_Enhancer
    chr1	1229469	1229579	NM_030649_exon_2_0_chr1_1229470_r	0	-	chr1	1228937	1229937	6_Weak_Enhancer
    chr1	1234724	1234736	NM_030649_exon_14_0_chr1_1234725_r	0	-	chr1	1234137	1234937	7_Weak_Enhancer
    chr1	1245060	1245231	NM_153339_exon_4_0_chr1_1245061_f	0	+	chr1	1244937	1245337	4_Strong_Enhancer

    bedtools intersect -a exons.bed \
                       -b <(grep Enhancer hesc.chromHmm.bed) \
                       -wa -wb -f 1.0 -u \
    | wc -l
    13746


5. What fraction of the GWAS SNPs are exonic?

Answer (Any idea why  we need -u?):

    wc -l gwas.bed
    17680 gwas.bed

    bedtools intersect -a gwas.bed -b exons.bed -u | wc -l
    1625

    echo "foo" | awk '{print 1625/17680}'
    0.0919118

6. What fraction of the GWAS SNPs are lie in either enhancers or promoters in the hESC data we have?

Answer (Any idea why  we need -u?):

    bedtools intersect -a gwas.bed -b <(egrep "Enhancer|Promoter" hesc.chromHmm.bed) -u \
    | wc -l
    1285

    echo "foo" | awk '{print 1285/17680}'
    0.072681

7. Create intervals representing the canonical 2bp splice sites on either side of each exon (don't worry about excluding splice sites at the first or last exon). (Hint - have a look at the [flank](http://bedtools.readthedocs.org/en/latest/content/tools/flank.html) tool.)

Answer:

    bedtools flank -l 2 -r 2 -i exons.bed -g genome.txt > splice-sites.bed

Or:

	bedtools slop -b 2 -i exons.bed -g genome.txt > exons.plus2.bed
	bedtools subtract -a exons.plus2.bed -b exons.bed > splice-sites.bed


8. What is the Jaccard statistic between CpG and hESC enhancers? Compare that to the Jaccard statistic between CpG and hESC promoters. Does the result make sense? (Hint - you will need `grep`).

Answer:

    bedtools jaccard -a cpg.bed -b <(grep Enhancer hesc.chromHmm.bed)
    intersection	union-intersection	jaccard	n_intersections
    1148180	132977386	0.0086344	4969

    bedtools jaccard -a cpg.bed -b <(grep Promoter hesc.chromHmm.bed)
    intersection	union-intersection	jaccard	n_intersections
    15661111	53551816	0.292448	20402


9. What would you expect the Jaccard statistic to look like if promoters were randomly distributed throughout the genome?  (Hint - you will need the [shuffle](http://bedtools.readthedocs.org/en/latest/content/tools/shuffle.html) tool.)

Answer:

    bedtools shuffle -i <(grep Promoter hesc.chromHmm.bed) -g genome.txt \
      | sort -k1,1 -k2,2n \
    > promoters.shuffled.bed

    bedtools jaccard -a cpg.bed -b promoters.shuffled.bed
    intersection	union-intersection	jaccard	n_intersections
    294071	68556207	0.00428949	78

10. Which hESC ChromHMM state (e.g., 11_Weak_Txn, 10_Txn_Elongation) represents the most number of base pairs in the genome? (Hint: you will need to use `awk` or `perl` here, as well as the [groupby](http://bedtools.readthedocs.org/en/latest/content/tools/groupby.html) tool.)
