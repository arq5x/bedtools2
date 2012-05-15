BT=../../bin/bedtools

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

rm *.bam