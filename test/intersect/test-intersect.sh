BT=${BT-../../bin/bedtools}

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

###########################################################
#  Test a basic self intersection
############################################################
echo "    intersect.t01...\c"
echo \
"chr1	10	20	a1	1	+
chr1	100	200	a2	2	-" > exp
$BT intersect -a a.bed -b a.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test a basic self intersection with -v
############################################################
echo "    intersect.t02...\c"
# expectation is empty set
touch exp
$BT intersect -a a.bed -b a.bed -v > obs
check obs exp
rm obs exp

###########################################################
#  Test -c
############################################################
echo "    intersect.t03...\c"
echo \
"chr1	10	20	a1	1	+	0
chr1	100	200	a2	2	-	2" > exp
$BT intersect -a a.bed -b b.bed -c > obs
check obs exp
rm obs exp

###########################################################
#  Test -c with -s
############################################################
echo "    intersect.t04...\c"
echo \
"chr1	10	20	a1	1	+	0
chr1	100	200	a2	2	-	1" > exp
$BT intersect -a a.bed -b b.bed -c -s > obs
check obs exp
rm obs exp

###########################################################
#  Test -c with -s and -f
############################################################
echo "    intersect.t05...\c"
echo \
"chr1	10	20	a1	1	+	0
chr1	100	200	a2	2	-	0" > exp
$BT intersect -a a.bed -b b.bed -c -s -f 0.1 > obs
check obs exp
rm obs exp

###########################################################
#  Test plain a and b intersect
############################################################
echo "    intersect.t06...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a.bed -b b.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test with -wa
############################################################
echo "    intersect.t07...\c"
echo \
"chr1	100	200	a2	2	-
chr1	100	200	a2	2	-" > exp
$BT intersect -a a.bed -b b.bed -wa > obs
check obs exp
rm obs exp

###########################################################
#  Test with -wa and -wb
############################################################
echo "    intersect.t08...\c"
echo \
"chr1	100	200	a2	2	-	chr1	90	101	b2	2	-
chr1	100	200	a2	2	-	chr1	100	110	b3	3	+" > exp
$BT intersect -a a.bed -b b.bed -wa -wb > obs
check obs exp
rm obs exp

###########################################################
#  Test with -wo (write overlap)
############################################################
echo "    intersect.t09...\c"
echo \
"chr1	100	200	a2	2	-	chr1	90	101	b2	2	-	1
chr1	100	200	a2	2	-	chr1	100	110	b3	3	+	10" > exp
$BT intersect -a a.bed -b b.bed -wo > obs
check obs exp
rm obs exp

###########################################################
#  Test with -wao (write all overlap)
############################################################
echo "    intersect.t10...\c"
echo \
"chr1	10	20	a1	1	+	.	-1	-1	.	-1	.	0
chr1	100	200	a2	2	-	chr1	90	101	b2	2	-	1
chr1	100	200	a2	2	-	chr1	100	110	b3	3	+	10" > exp
$BT intersect -a a.bed -b b.bed -wao > obs
check obs exp
rm obs exp

###########################################################
#  Test with -wo (write overlap) with -s
############################################################
echo "    intersect.t11...\c"
echo \
"chr1	100	200	a2	2	-	chr1	90	101	b2	2	-	1" > exp
$BT intersect -a a.bed -b b.bed -wo -s > obs
check obs exp
rm obs exp


###########################################################
#  Test with -wao (write all overlap) with -s
############################################################
echo "    intersect.t12...\c"
echo \
"chr1	10	20	a1	1	+	.	-1	-1	.	-1	.	0
chr1	100	200	a2	2	-	chr1	90	101	b2	2	-	1" > exp
$BT intersect -a a.bed -b b.bed -wao -s > obs
check obs exp
rm obs exp

###########################################################
#  Test A as -
############################################################
echo "    intersect.t13...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a.bed | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test A as stdin
############################################################
echo "    intersect.t14...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a.bed | $BT intersect -a stdin -b b.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test B as -
############################################################
echo "    intersect.t15...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat b.bed | $BT intersect -a a.bed -b - > obs
check obs exp
rm obs exp

###########################################################
#  Test A as stdin
############################################################
echo "    intersect.t16...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat b.bed | $BT intersect -a a.bed -b stdin > obs
check obs exp
rm obs exp


###########################################################
###########################################################
#                       -split                            #
###########################################################
###########################################################
samtools view -Sb one_block.sam > one_block.bam 2>/dev/null
samtools view -Sb two_blocks.sam > two_blocks.bam 2>/dev/null
samtools view -Sb three_blocks.sam > three_blocks.bam 2>/dev/null


##################################################################
#  Test three blocks matches BED without -split
##################################################################
echo "    intersect.t17...\c"
echo \
"three_blocks	16	chr1	1	40	10M10N10M10N10M	*	0	0	GAAGGCCACCGCCGCGGTTATTTTCCTTCA	CCCDDB?=FJIIJIGFJIJHIJJJJJJJJI	MD:Z:50" > exp
$BT intersect -abam three_blocks.bam -b three_blocks_nomatch.bed | samtools view - > obs
check obs exp
rm obs exp

##################################################################
#  Test three blocks does not match BED with -split
##################################################################
echo "    intersect.t18...\c"
touch exp
$BT intersect -abam three_blocks.bam -b three_blocks_nomatch.bed -split | samtools view - > obs
check obs exp
rm obs exp

##################################################################
#  Test three blocks matches BED with -split
##################################################################
echo "    intersect.t19...\c"
echo \
"three_blocks	16	chr1	1	40	10M10N10M10N10M	*	0	0	GAAGGCCACCGCCGCGGTTATTTTCCTTCA	CCCDDB?=FJIIJIGFJIJHIJJJJJJJJI	MD:Z:50" > exp
$BT intersect -abam three_blocks.bam -b three_blocks_match.bed -split | samtools view - > obs
check obs exp
rm obs exp

##################################################################
#  Test three blocks does not match BED with -split and -s
#  BAM os -, BED is +
##################################################################
echo "    intersect.t20...\c"
touch exp
$BT intersect -abam three_blocks.bam -b three_blocks_match.bed -split -s | samtools view - > obs
check obs exp
rm obs exp

##################################################################
#  Test three blocks  match BED that overlap 1bp with -split
##################################################################
echo "    intersect.t21...\c"
echo \
"three_blocks	16	chr1	1	40	10M10N10M10N10M	*	0	0	GAAGGCCACCGCCGCGGTTATTTTCCTTCA	CCCDDB?=FJIIJIGFJIJHIJJJJJJJJI	MD:Z:50" > exp
$BT intersect -abam three_blocks.bam -b three_blocks_match_1bp.bed -split | samtools view - > obs
check obs exp
rm obs exp

################################################################################
#  Test three blocks does not match BED that overlap 1bp with -split and -f 0.1
################################################################################
echo "    intersect.t22...\c"
touch exp
$BT intersect -abam three_blocks.bam -b three_blocks_match_1bp.bed -split -f 0.1 | samtools view - > obs
check obs exp
rm obs exp



###########################################################
#  Test three blocks with -split -wo only shows 5 overlap
#  bases, not ten.
############################################################
echo "    intersect.t22.a...\c"
echo "chr1	0	50	three_blocks_match	0	+	0	0	0	3	10,10,10,	0,20,40,	chr1	5	15	5" > exp
$BT intersect -a three_blocks_match.bed -b d.bed -split -wo > obs
check obs exp
rm obs exp

###########################################################
#  Same test but for BAM file
############################################################
echo "    intersect.t22.b...\c"
echo "chr1	0	50	three_blocks_match	255	+	0	50	0,0,0	3	10,10,10,	0,20,40,	chr1	5	15	5" > exp
$BT intersect -a three_blocks_match.bam -b d.bed -split -wo -bed > obs
check obs exp
rm obs exp


###########################################################
#  Test three blocks with -split -wo, and DB record also has
#  blocks that somewhat intersect
############################################################
echo "    intersect.t22.c...\c"
echo "chr1	0	50	three_blocks_match	0	+	0	0	0	3	10,10,10,	0,20,40,	chr1	0	45	three_blocks_match	0	+	0	0	0	2	5,10,	25,35,	10" > exp
$BT intersect -a three_blocks_match.bed -b two_blocks_partial.bed -split -wo > obs
check obs exp
rm obs exp

###########################################################
#  Same test but for BAM file
############################################################
echo "    intersect.t22.d...\c"
echo "chr1	0	50	three_blocks_match	255	+	0	50	0,0,0	3	10,10,10,	0,20,40,	chr1	0	45	three_blocks_match	0	+	0	0	0	2	5,10,	25,35,	10" > exp
$BT intersect -a three_blocks_match.bam -b two_blocks_partial.bed -split -wo -bed > obs
check obs exp
rm obs exp

###########################################################
#  Test three blocks with -split -wo, and DB record also has
#  blocks that do not intersect
############################################################
echo "    intersect.t22.e...\c"
touch exp
$BT intersect -a three_blocks_match.bed -b three_blocks_nomatch.bed -split -wo > obs
check obs exp
rm obs exp


###########################################################
#  Same test but for BAM file
############################################################
echo "    intersect.t22.f...\c"
touch exp
$BT intersect -a three_blocks_match.bam -b three_blocks_nomatch.bed -split -wo -bed > obs
check obs exp
rm obs exp



###########################################################
#  Added split test from bug150. See that overlap bases
# with -split are correctly affected by overlap fraction
############################################################
echo "    intersect.t22.g...\c"
echo \
"chr2	1000	16385	A	0	-	0	0	0	2	1,1,	0,15384,	chr2	1000	16385	A	0	-	0	0	0	2	1,1,	0,15384,	2" > exp
$BT intersect -a bug150_a.bed -b bug150_b.bed -s -split -wo > obs
check exp obs
rm exp obs


##################################################################
#  Test that only the mapped read is is found as an intersection
##################################################################
echo "    intersect.t23...\c"
echo \
"mapped	16	chr1	1	40	30M	*	0	0	GAAGGCCACCGCCGCGGTTATTTTCCTTCA	CCCDDB?=FJIIJIGFJIJHIJJJJJJJJI	MD:Z:50" > exp
samtools view -Sb mapped_and_unmapped.sam 2>/dev/null | $BT intersect -abam - -b a.bed | samtools view - > obs
check obs exp
rm obs exp

##################################################################
#  Test that an unmapped read is handled properly with -v
##################################################################
echo "    intersect.t24...\c"
echo \
"umapped	4	*	1	40	30M	*	0	0	GAAGGCCACCGCCGCGGTTATTTTCCTTCA	CCCDDB?=FJIIJIGFJIJHIJJJJJJJJI	MD:Z:50" > exp
samtools view -Sb mapped_and_unmapped.sam 2>/dev/null | $BT intersect -abam - -b a.bed -v | samtools view - > obs
check obs exp
rm obs exp

##################################################################
#  Test -c with BAM input
##################################################################
echo "    intersect.t25...\c"
echo \
"chr1	0	30	one_blocks	40	-	0	30	0,0,0	1	30,	0,	1" > exp
$BT intersect -abam one_block.bam -b c.bed -c -bed > obs
check obs exp
rm obs exp

##################################################################
#  Test -wo with BAM input
##################################################################
echo "    intersect.t26...\c"
echo \
"chr1	0	30	one_blocks	40	-	0	30	0,0,0	1	30,	0,	chr1	0	100	c1	1	+	30" > exp
$BT intersect -abam one_block.bam -b c.bed -wo -bed > obs
check obs exp
rm obs exp

##################################################################
#  Test -wao with BAM input
##################################################################
echo "    intersect.t27...\c"
echo \
"chr1	0	30	one_blocks	40	-	0	30	0,0,0	1	30,	0,	chr1	0	100	c1	1	+	30" > exp
$BT intersect -abam one_block.bam -b c.bed -wo -bed > obs
check obs exp

##################################################################
#  Test BED3 with BED3
##################################################################
echo "    intersect.t28...\c"
echo \
"chr1^I10^I20^Ichr1^I10^I20" > exp
$BT intersect -a bed3.bed -b bed3.bed -wa -wb | cat -t > obs
check obs exp

##################################################################
#  Test BED4 with BED3
##################################################################
echo "    intersect.t29...\c"
echo \
"chr1^I10^I20^I345.7^Ichr1^I10^I20" > exp
$BT intersect -a bed4.bed -b bed3.bed -wa -wb | cat -t > obs
check obs exp

##################################################################
#  Test BED5 with BED3
##################################################################
echo "    intersect.t30...\c"
echo \
"chr1^I10^I20^I345.7^Iwhy?^Ichr1^I10^I20" > exp
$BT intersect -a bed5.bed -b bed3.bed -wa -wb | cat -t > obs
check obs exp

##################################################################
#  Test BED6 (without a proper strand) with BED3
##################################################################
echo "    intersect.t31...\c"
echo \
"chr1^I10^I20^I345.7^Iwhy?^I11^Ichr1^I10^I20" > exp
$BT intersect -a bed6.bed -b bed3.bed -wa -wb | cat -t > obs
check obs exp

##################################################################
#  Test BED6 (with a strand) with BED3
##################################################################
echo "    intersect.t32...\c"
echo \
"chr1^I10^I20^I345.7^Iwhy?^I-^Ichr1^I10^I20" > exp
$BT intersect -a bed6.strand.bed -b bed3.bed -wa -wb | cat -t > obs
check obs exp

##################################################################
#  Test BED PLUS with BED3
##################################################################
echo "    intersect.t33...\c"
echo \
"chr1^I10^I20^I345.7^Iwhy?^I11^Ifoo^Ibar^Ibiz^I11^Ibang^Ibop^I99^Ichr1^I10^I20" > exp
$BT intersect -a bedplus.bed -b bed3.bed -wa -wb | cat -t > obs
check obs exp

##################################################################
#  Test for strand matches with BED3
##################################################################
echo "    intersect.t34...\c"
echo \
"chr1^I10^I20^I345.7^Iwhy?^I11^Ichr1^I10^I20^I345.7^Iwhy?^I-" > exp
$BT intersect -a bed6.bed -b bed6.strand.bed -wa -wb | cat -t > obs
check obs exp

##################################################################
#  Test for strand matches with BED3
##################################################################
echo "    intersect.t35...\c"
echo \
"chr1^I10^I20^I345.7^Iwhy?^I-^Ichr1^I10^I20^I345.7^Iwhy?^I-" > exp
$BT intersect -a bed6.strand.bed -b bed6.strand2.bed -wa -wb -s | cat -t > obs
check obs exp

##################################################################
#  Test for strand matches with BED3
##################################################################
echo "    intersect.t36...\c"
echo \
"chr1^I10^I20^I345.7^Iwhy?^I-^Ichr1^I11^I21^I345.7^Iwhy?^I+" > exp
$BT intersect -a bed6.strand.bed -b bed6.strand2.bed -wa -wb -S | cat -t > obs
check obs exp
rm obs exp

##################################################################
#  Test that intersect of bed query with BAM DB gives Bed output.
##################################################################
echo "    intersect.t37...\c"
echo \
"chr1	10	20	a1	1	+
chr1	100	200	a2	2	-" > exp
$BT intersect -a a.bed -b a.bam > obs
check obs exp
rm obs exp

##################################################################
#  Test that -split works on identical records even if
# -f 1 is ued (100% overlap)
##################################################################
echo "    intersect.t38...\c"
echo \
"chr1	1	100	A	0	+	1	100	0	2	10,10	0,90	chr1	1	100	B	0	+	1	100	0	2	10,10	0,90	20" > exp
$BT intersect -wao  -f 1 -split -a splitBug155_a.bed -b splitBug155_b.bed > obs
check obs exp
rm obs exp

##################################################################
#  Test that fractional overlap must be greater than 0.0
##################################################################
echo "    intersect.t39...\c"
echo \
"***** ERROR: -f must be in the range (0.0, 1.0]. *****" > exp
$BT intersect -a a.bed -b b.bed -f 0.0 2>&1 > /dev/null | cat - | head -2 | tail -1 > obs
check exp obs
rm exp obs

##################################################################
#  Test that fractional overlap must be <= than 1.0
##################################################################
echo "    intersect.t40...\c"
echo \
"***** ERROR: -f must be in the range (0.0, 1.0]. *****" > exp
$BT intersect -a a.bed -b b.bed -f 1.00001 2>&1 > /dev/null | cat - | head -2 | tail -1 > obs
check exp obs
rm exp obs

##################################################################
#  bug167_strandSweep.bed
##################################################################
echo "    intersect.t41...\c"
echo \
"22" > exp
$BT intersect -a bug167_strandSweep.bed -b bug167_strandSweep.bed -sorted -s -wa -wb | wc -l > obs
sed -i 's/^\s*//' obs
check exp obs
rm exp obs

##################################################################
#  bug167_strandSweep.bed
##################################################################
echo "    intersect.t42...\c"
echo \
"20" > exp
$BT intersect -a bug167_strandSweep.bed -b bug167_strandSweep.bed -sorted -S -wa -wb | wc -l > obs
sed -i 's/^\s*//' obs
check exp obs
rm exp obs

rm one_block.bam two_blocks.bam three_blocks.bam


##################################################################
# Bug 187 0 length records
##################################################################
echo "    intersect.t43...\c"
echo \
"chr7	33059403	33059403	chr7	33059336	33060883	NT5C3A	intron	0
chr7	33059403	33059403	chr7	32599076	33069221	NAq	intron	0" > exp
$BT intersect -a bug187_a.bed -b bug187_b.bed -wao > obs
check exp obs
rm exp obs

##################################################################
# see that naming conventions are tested with unsorted data.
##################################################################
echo "    intersect.t44...\c"
echo \
"***** WARNING: File nonamecheck_a.bed has a record where naming convention (leading zero) is inconsistent with other files:
chr1	10	20" > exp
$BT intersect -a nonamecheck_a.bed -b nonamecheck_b.bed 2>&1 > /dev/null | cat - | head -2 > obs
check exp obs
rm exp obs


##################################################################
# see that differently named chroms don't work with -sorted
##################################################################
echo "    intersect.t45...\c"
echo \
"***** WARNING: File nonamecheck_b.bed has a record where naming convention (leading zero) is inconsistent with other files:
chr01	15	25" > exp
$BT intersect -a nonamecheck_a.bed -b nonamecheck_b.bed -sorted 2>&1 > /dev/null | cat - | head -2 > obs
check exp obs
rm exp obs

##################################################################
# see that differently named chroms  -sorted and -nonamecheck
# don't complain with -nonamecheck
##################################################################
echo "    intersect.t46...\c"
touch exp
$BT intersect -a nonamecheck_a.bed -b nonamecheck_b.bed -sorted -nonamecheck 2>&1 > /dev/null | cat - > obs
check exp obs
rm exp obs

##################################################################
# see that SVLEN in VCF files is treated as zero length
# records when the SV type is an insertion
##################################################################
echo "    intersect.t47...\c"
echo \
"chr1	1	a	G	<DEL>	70.90
chr1	1	a	G	<DEL>	70.90
chr1	4	a	G	<INS>	70.90" > exp
$BT intersect -a bug223_sv1_a.vcf -b bug223_sv1_b.vcf | cut -f1-6 > obs
check exp obs
rm exp obs

##################################################################
# see that SVLEN in VCF files can handle multiple numbers,
# at end of line, followed by NULL.
##################################################################
echo "    intersect.t48...\c"
echo \
"chr1	1	a	G	<DEL>	70.90
chr1	1	a	G	<DEL>	70.90
chr1	4	a	G	<DEL>	70.90
chr1	4	a	G	<DEL>	70.90" > exp
$BT intersect -a bug223_d.vcf -b bug223_d.vcf | cut -f1-6 > obs
check exp obs
rm exp obs

##################################################################
# see that SVLEN in VCF files can handle multiple numbers,
# at end of line, followed by a tab
##################################################################
echo "    intersect.t49...\c"
echo \
"chr1	1	a	G	<DEL>	70.90
chr1	1	a	G	<DEL>	70.90
chr1	4	a	G	<DEL>	70.90
chr1	4	a	G	<DEL>	70.90" > exp
$BT intersect -a bug223_e.vcf -b bug223_e.vcf | cut -f1-6 > obs
check exp obs
rm exp obs

##################################################################
# see that SVLEN in VCF files can handle single numbers,
# at end of line, followed by null
##################################################################
echo "    intersect.t50...\c"
echo \
"chr1	1	a	G	<DEL>	70.90
chr1	1	a	G	<DEL>	70.90
chr1	4	a	G	<DEL>	70.90
chr1	4	a	G	<DEL>	70.90" > exp
$BT intersect -a bug223_f.vcf -b bug223_f.vcf | cut -f1-6 > obs
check exp obs
rm exp obs

##################################################################
# Bug 44: test that bgzipped vcf file works correctly
# with race condition
##################################################################
echo "    intersect.t51...\c"
echo \
"MT	2706	.	A	G	2965	PASS	BRF=0.05;FR=1;HP=1;HapScore=1;MGOF=17;MMLQ=30;MQ=62.05;NF=7607;NR=8147;PP=2965;QD=20;SC=AGGCGGGCATAACACAGCAAG;SbPval=0.52;Source=Platypus;TC=15840;TCF=7679;TCR=8161;TR=15754;WE=2749;WS=2693;CSQ=G|ENSG00000198763|ENST00000361453|Transcript|upstream_gene_variant||||||rs2854128|1764|1|MT-ND2|HGNC|7456|protein_coding|YES||ENSP00000355046|NU2M_HUMAN|Q7GXY9_HUMAN&Q5Q3P5_HUMAN&Q14X33_HUMAN&Q14WT3_HUMAN&A6ZH82_HUMAN&A6ZGN8_HUMAN&A6ZGG3_HUMAN|UPI0000000AA2||||||A:0.1656|||||||||||||,G|ENSG00000210151|ENST00000387416|Transcript|downstream_gene_variant||||||rs2854128|4740|-1|MT-TS1|HGNC|7497|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210077|ENST00000387342|Transcript|downstream_gene_variant||||||rs2854128|1036|1|MT-TV|HGNC|7500|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210144|ENST00000387409|Transcript|downstream_gene_variant||||||rs2854128|3120|-1|MT-TY|HGNC|7502|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210117|ENST00000387382|Transcript|upstream_gene_variant||||||rs2854128|2806|1|MT-TW|HGNC|7501|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210107|ENST00000387372|Transcript|downstream_gene_variant||||||rs2854128|1623|-1|MT-TQ|HGNC|7495|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210140|ENST00000387405|Transcript|downstream_gene_variant||||||rs2854128|3055|-1|MT-TC|HGNC|7477|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000211459|ENST00000389680|Transcript|downstream_gene_variant||||||rs2854128|1105|1|MT-RNR1|HGNC|7470|Mt_rRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210082|ENST00000387347|Transcript|non_coding_transcript_exon_variant&non_coding_transcript_variant|1036|||||rs2854128||1|MT-RNR2|HGNC|7471|Mt_rRNA|YES||||||||1/1|||A:0.1656|||||||||||||,G|ENSG00000210127|ENST00000387392|Transcript|downstream_gene_variant||||||rs2854128|2881|-1|MT-TA|HGNC|7475|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000198712|ENST00000361739|Transcript|upstream_gene_variant||||||rs2854128|4880|1|MT-CO2|HGNC|7421|protein_coding|YES||ENSP00000354876|COX2_HUMAN|Q7GXZ8_HUMAN&Q4R1L5_HUMAN&Q4R1L3_HUMAN&Q14XT3_HUMAN&K7WVJ5_HUMAN&H9E7W2_HUMAN&H9E7T7_HUMAN&H9E7P8_HUMAN&H9E7F7_HUMAN&E2DTL8_HUMAN&D3WYY9_HUMAN&D2Y6Y2_HUMAN&D2Y6Y1_HUMAN&B2YKU2_HUMAN|UPI0000000AA4||||||A:0.1656|||||||||||||,G|ENSG00000210049|ENST00000387314|Transcript|downstream_gene_variant||||||rs2854128|2059|1|MT-TF|HGNC|7481|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000198888|ENST00000361390|Transcript|upstream_gene_variant||||||rs2854128|601|1|MT-ND1|HGNC|7455|protein_coding|YES||ENSP00000354687|NU1M_HUMAN|Q85KV6_HUMAN&Q8WCX9_HUMAN&Q5Q757_HUMAN&Q14WI3_HUMAN&G3EBI1_HUMAN&D2Y6X8_HUMAN&D2Y6X6_HUMAN&A6ZHG8_HUMAN|UPI0000000AA1||||||A:0.1656|||||||||||||,G|ENSG00000209082|ENST00000386347|Transcript|upstream_gene_variant||||||rs2854128|524|1|MT-TL1|HGNC|7490|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000198804|ENST00000361624|Transcript|upstream_gene_variant||||||rs2854128|3198|1|MT-CO1|HGNC|7419|protein_coding|YES||ENSP00000354499|COX1_HUMAN|Q957U9_HUMAN&Q7GXY8_HUMAN&M9Z2G2_HUMAN&Q8HBX8_HUMAN&Q5Q1W2_HUMAN&Q4R1L4_HUMAN&Q14XD3_HUMAN&Q14X83_HUMAN&F8U4W0_HUMAN&D3WYY6_HUMAN&D3WYY5_HUMAN&D3WYY4_HUMAN&D2Y6W4_HUMAN&C8YAE4_HUMAN&C3UPN2_HUMAN&B7TCT8_HUMAN&B2Y9D8_HUMAN&A5YMT3_HUMAN&A1XP63_HUMAN&A0S1I7_HUMAN|UPI0000000AA3||||||A:0.1656|||||||||||||,G|ENSG00000210154|ENST00000387419|Transcript|upstream_gene_variant||||||rs2854128|4812|1|MT-TD|HGNC|7478|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210112|ENST00000387377|Transcript|upstream_gene_variant||||||rs2854128|1696|1|MT-TM|HGNC|7492|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210135|ENST00000387400|Transcript|downstream_gene_variant||||||rs2854128|2951|-1|MT-TN|HGNC|7493|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||,G|ENSG00000210100|ENST00000387365|Transcript|upstream_gene_variant||||||rs2854128|1557|1|MT-TI|HGNC|7488|Mt_tRNA|YES|||||||||||A:0.1656|||||||||||||;GR=3.07;PH=0.654;PS=0.002	GT:GL:GOF:GQ:NR:NV	1/1:-300,-298.01,0:3:99:2733:2718	1/1:-300,-298.01,0:17:99:6509:6461	1/1:-300,-298.01,0:2:99:6598:6575	MT	2591	2747	rRNA" > exp
$BT intersect -a bug44_a.vcf.gz -b bug44_b.bed -wa -wb > obs
check exp obs
rm exp obs

##################################################################
# Test basic -f functionality
##################################################################
echo "    intersect.t52...\c"
echo "chr1	10	12	a1	1	+
chr2	10	12	a2	1	-" > exp
$BT intersect -a x.bed -b y.bed -f 0.2 > obs
check exp obs
rm exp obs

echo "    intersect.t53...\c"
echo -n "" > exp
$BT intersect -a x.bed -b y.bed -f 0.21 > obs
check exp obs
rm exp obs

##################################################################
# Test basic -F functionality
##################################################################
echo "    intersect.t54...\c"
echo "chr1	10	20	a1	1	+	chr1	8	12	b1	1	+
chr2	10	20	a2	1	-	chr2	8	12	b2	1	+" > exp
$BT intersect -a x.bed -b y.bed -F 0.21 -wa -wb > obs
check exp obs
rm exp obs

##################################################################
# Test basic -f with -F
##################################################################
echo "    intersect.t55...\c"
echo -n "" > exp
$BT intersect -a x.bed -b y.bed -f 0.21 -F 0.21 -wa -wb > obs
check exp obs
rm exp obs

echo "    intersect.t56...\c"
echo -n "" > exp
$BT intersect -a x.bed -b y.bed -f 0.21 -r -wa -wb > obs
check exp obs
rm exp obs

echo "    intersect.t57...\c"
echo "chr1	10	20	a1	1	+	chr1	8	12	b1	1	+
chr2	10	20	a2	1	-	chr2	8	12	b2	1	+" > exp
$BT intersect -a x.bed -b y.bed -f 0.19 -r -wa -wb > obs
check exp obs
rm exp obs

echo "    intersect.t58...\c"
echo "chr1	10	20	a1	1	+	chr1	8	12	b1	1	+
chr2	10	20	a2	1	-	chr2	8	12	b2	1	+" > exp
$BT intersect -a x.bed -b y.bed -F 0.50 -wa -wb > obs
check exp obs
rm exp obs

echo "    intersect.t59...\c"
echo "chr1	10	20	a1	1	+	chr1	8	12	b1	1	+
chr2	10	20	a2	1	-	chr2	8	12	b2	1	+" > exp
$BT intersect -a x.bed -b y.bed -f 0.19 -F 0.21 -wa -wb > obs
check exp obs
rm exp obs

echo "    intersect.t60...\c"
echo "chr1	10	20	a1	1	+	chr1	8	12	b1	1	+
chr2	10	20	a2	1	-	chr2	8	12	b2	1	+" > exp
$BT intersect -a x.bed -b y.bed -f 0.19 -F 0.21 -wa -wb > obs
check exp obs
rm exp obs

echo "    intersect.t61...\c"
echo "chr1	10	20	a1	1	+	chr1	8	12	b1	1	+
chr2	10	20	a2	1	-	chr2	8	12	b2	1	+" > exp
$BT intersect -a x.bed -b y.bed -f 0.19 -F 0.50 -wa -wb > obs
check exp obs
rm exp obs

echo "    intersect.t62...\c"
echo -n "" > exp
$BT intersect -a x.bed -b y.bed -f 0.19 -F 0.51 -wa -wb > obs
check exp obs
rm exp obs

##################################################################
# Test basic -f with -F with BAM
##################################################################
echo "    intersect.t63...\c"
echo -n "" > exp
$BT intersect -a x.bam -b y.bed -f 0.21 -F 0.21 -wa -wb | samtools view - > obs
check exp obs
rm exp obs

echo "    intersect.t64...\c"
echo "a1	0	chr1	11	255	10M	*	0	0	*	*
a2	16	chr2	11	255	10M	*	0	0	*	*" > exp
$BT intersect -a x.bam -b y.bed -f 0.19 -F 0.21 -wa -wb | samtools view - > obs
check exp obs
rm exp obs

echo "    intersect.t65...\c"
echo "a1	0	chr1	11	255	10M	*	0	0	*	*
a2	16	chr2	11	255	10M	*	0	0	*	*" > exp
$BT intersect -a x.bam -b y.bed -f 0.19 -F 0.21 -wa -wb | samtools view - > obs
check exp obs
rm exp obs

echo "    intersect.t66...\c"
echo "a1	0	chr1	11	255	10M	*	0	0	*	*
a2	16	chr2	11	255	10M	*	0	0	*	*" > exp
$BT intersect -a x.bam -b y.bed -f 0.19 -F 0.50 -wa -wb | samtools view - > obs
check exp obs
rm exp obs

echo "    intersect.t67...\c"
echo -n "" > exp
$BT intersect -a x.bam -b y.bed -f 0.19 -F 0.51 -wa -wb | samtools view - > obs
check exp obs
rm exp obs


##################################################################
# Test basic -f with -F and -e
##################################################################
echo "    intersect.t68...\c"
echo "chr1	10	20	a1	1	+	chr1	8	12	b1	1	+
chr2	10	20	a2	1	-	chr2	8	12	b2	1	+" > exp
$BT intersect -a x.bed -b y.bed -f 0.21 -F 0.21 -wa -wb -e > obs
check exp obs
rm exp obs


##################################################################
# Test basic -split with BED12 w/ and w/o trailing commas for blocks
# Issue 366
##################################################################
echo "    intersect.t69...\c"
echo "chr1	0	45	oneblock_comma	0	+	0	0	0	1	45,	0,
chr1	0	45	oneblock_comma	0	+	0	0	0	1	45,	0,
chr1	0	45	oneblock_comma	0	+	0	0	0	1	45,	0,
chr1	0	45	oneblock_comma	0	+	0	0	0	1	45,	0,
chr1	0	45	oneblock_nocomma	0	+	0	0	0	1	45	0
chr1	0	45	oneblock_nocomma	0	+	0	0	0	1	45	0
chr1	0	45	oneblock_nocomma	0	+	0	0	0	1	45	0
chr1	0	45	oneblock_nocomma	0	+	0	0	0	1	45	0
chr1	0	45	three_blocks_comma	0	+	0	0	0	3	10,10,10,	0,20,40,
chr1	0	45	three_blocks_comma	0	+	0	0	0	3	10,10,10,	0,20,40,
chr1	0	50	three_blocks_comma	0	+	0	0	0	3	10,10,10,	0,20,40,
chr1	0	50	three_blocks_comma	0	+	0	0	0	3	10,10,10,	0,20,40,
chr1	0	45	three_blocks_nocomma	0	+	0	0	0	3	10,10,10	0,20,40
chr1	0	45	three_blocks_nocomma	0	+	0	0	0	3	10,10,10	0,20,40
chr1	0	50	three_blocks_nocomma	0	+	0	0	0	3	10,10,10	0,20,40
chr1	0	50	three_blocks_nocomma	0	+	0	0	0	3	10,10,10	0,20,40" > exp
$BT intersect  -a blocks.bed12 -b blocks.bed12 -split > obs
check exp obs
rm exp obs

##################################################################
# Issue 311
##################################################################
echo "    intersect.t70...\c"
echo "1	31	32	1	32	.	A	T	0	PASS	DP=22" > exp
$BT intersect -a a.issue311.bed -b b.issue311.vcf -wb > obs
check exp obs
rm exp obs

echo "    intersect.t71...\c"
echo "1	31	32	1	pseudogene	exon	32	32	.	-	.	 gene_id \"ENSG00000224777\"; transcript_id \"ENST00000424047\"; exon_number \"1\"; gene_name \"OR4F2P\"; transcript_name \"OR4F2P-001\";" > exp
$BT intersect -a a.issue311.bed -b b.issue311.gff -wb > obs
check exp obs
rm exp obs

# GFF   ----------
# BED        -----------
echo "    intersect.t72...\c"
echo "1	.	.	16	20	.	-	.	." > exp
$BT intersect -a <(echo -e  "1\t.\t.\t10\t20\t.\t-\t.\t.") -b <(echo -e  "1\t15\t25") > obs
check exp obs
rm exp obs

# BED   ----------
# GFF        -----------
echo "    intersect.t73...\c"
echo "1	15	20" > exp
$BT intersect -a <(echo -e  "1\t15\t25") -b <(echo -e  "1\t.\t.\t10\t20\t.\t-\t.\t.")  > obs
check exp obs
rm exp obs

# GFF        ----------
# BED    -----------
echo "    intersect.t74...\c"
echo "1	.	.	15	20	.	-	.	." > exp
$BT intersect -a <(echo -e  "1\t.\t.\t15\t25\t.\t-\t.\t.") -b <(echo -e  "1\t10\t20") > obs
check exp obs
rm exp obs

# BED        ----------
# GFF    -----------
echo "    intersect.t75...\c"
echo "1	15	20" > exp
$BT intersect  -a <(echo -e  "1\t15\t25") -b <(echo -e  "1\t.\t.\t10\t20\t.\t-\t.\t.") > obs
check exp obs
rm exp obs

# GFF        ----------
# BED    -------------------
echo "    intersect.t76...\c"
echo "1	.	.	15	20	.	-	.	." > exp
$BT intersect -a <(echo -e  "1\t.\t.\t15\t20\t.\t-\t.\t.") -b <(echo -e  "1\t10\t25") > obs
check exp obs
rm exp obs

# BED        ----------
# GFF    --------------------
echo "    intersect.t77...\c"
echo "1	15	20" > exp
$BT intersect -a <(echo -e "1\t15\t20") -b <(echo -e "1\t.\t.\t10\t25\t.\t-\t.\t.")  > obs
check exp obs
rm exp obs

# GFF        ----------
# BED    -------------------
echo "    intersect.t78...\c"
echo "1	.	.	16	20	.	-	.	." > exp
$BT intersect -a <(echo -e "1\t.\t.\t10\t25\t.\t-\t.\t.") -b <(echo -e "1\t15\t20") > obs
check exp obs
rm exp obs

# BED        ----------
# GFF    --------------------
echo "    intersect.t79...\c"
echo "1	14	20" > exp
$BT intersect -a <(echo -e "1\t10\t25") -b <(echo -e "1\t.\t.\t15\t20\t.\t-\t.\t.")  > obs
check exp obs
rm exp obs


##################################################################
# Issue 316. SVLEN and END
##################################################################
echo "    intersect.t72...\c"

echo "##fileformat=VCF4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	32	.	A	<DEL>	0	PASS	DP=22;END=52" > a
echo "##fileformat=VCF4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	52	.	A	T	0	PASS	DP=22;SVLEN=100" > b
echo "1	32	.	A	<DEL>	0	PASS	DP=22;END=52" > exp
$BT intersect -a a -b b > obs
check exp obs
rm exp obs a b

echo "    intersect.t73...\c"

echo "##fileformat=VCF4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	32	.	A	<DEL>	0	PASS	DP=22;END=51" > a
echo "##fileformat=VCF4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	52	.	A	T	0	PASS	DP=22;SVLEN=100" > b
echo -n "" > exp
$BT intersect -a a -b b > obs
check exp obs
rm exp obs a b



cd multi_intersect
bash test-multi_intersect.sh
cd ..

cd sortAndNaming
bash test-sort-and-naming.sh
cd ..
