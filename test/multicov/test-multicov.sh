set -e;
BT=${BT-../../bin/bedtools}
htsutil=${htsutil-../htsutil}

FAILURES=0;

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	FAILURES=$(expr $FAILURES + 1);
		echo fail
	fi
}

$htsutil samtoindexedbam one_block.sam one_block.bam
$htsutil samtoindexedbam two_blocks.sam two_blocks.bam
$htsutil samtoindexedbam test-multi.sam test-multi.bam
$htsutil samtoindexedbam test-multi.2.sam test-multi.2.bam

##################################################################
#  Test one block matches all BEDs
##################################################################
echo -e "    multicov.t1...\c"
echo \
"chr1	15	20	a1	1	+	1
chr1	15	27	a2	2	+	1
chr1	15	20	a3	3	-	1
chr1	15	27	a4	4	-	1" > exp
$BT multicov -bams one_block.bam -bed multicov.bed > obs
check obs exp
rm obs exp

##################################################################
#  Test one block matches based on _same_ strand
##################################################################
echo -e "    multicov.t2...\c"
echo \
"chr1	15	20	a1	1	+	0
chr1	15	27	a2	2	+	0
chr1	15	20	a3	3	-	1
chr1	15	27	a4	4	-	1" > exp
$BT multicov -bams one_block.bam -bed multicov.bed -s > obs
check obs exp
rm obs exp

##################################################################
#  Test one block matches based on _different_ strands
##################################################################
echo -e "    multicov.t3...\c"
echo \
"chr1	15	20	a1	1	+	1
chr1	15	27	a2	2	+	1
chr1	15	20	a3	3	-	0
chr1	15	27	a4	4	-	0" > exp
$BT multicov -bams one_block.bam -bed multicov.bed -S > obs
check obs exp
rm obs exp

##################################################################
#  Test split alignment matches 
##################################################################
echo -e "    multicov.t4...\c"
echo \
"chr1	15	20	a1	1	+	1
chr1	15	27	a2	2	+	1
chr1	15	20	a3	3	-	1
chr1	15	27	a4	4	-	1" > exp
$BT multicov -bams two_blocks.bam -bed multicov.bed > obs
check obs exp
rm obs exp

##################################################################
#  Test split alignment matches with -split
##################################################################
echo -e "    multicov.t5...\c"
echo \
"chr1	15	20	a1	1	+	0
chr1	15	27	a2	2	+	1
chr1	15	20	a3	3	-	0
chr1	15	27	a4	4	-	1" > exp
$BT multicov -bams two_blocks.bam -bed multicov.bed -split > obs
check obs exp
rm obs exp

##################################################################
#  Test split alignment matches with -split and -s
##################################################################
echo -e "    multicov.t6...\c"
echo \
"chr1	15	20	a1	1	+	0
chr1	15	27	a2	2	+	0
chr1	15	20	a3	3	-	0
chr1	15	27	a4	4	-	1" > exp
$BT multicov -bams two_blocks.bam -bed multicov.bed -split -s > obs
check obs exp
rm obs exp

##################################################################
#  Test split alignment matches with -split and -S
##################################################################
echo -e "    multicov.t7...\c"
echo \
"chr1	15	20	a1	1	+	0
chr1	15	27	a2	2	+	1
chr1	15	20	a3	3	-	0
chr1	15	27	a4	4	-	0" > exp
$BT multicov -bams two_blocks.bam -bed multicov.bed -split -S > obs
check obs exp
rm obs exp

##################################################################
#  Test split alignment matches with -split and -f
##################################################################
echo -e "    multicov.t8...\c"
echo \
"chr1	15	20	a1	1	+	0
chr1	15	27	a2	2	+	1
chr1	15	20	a3	3	-	0
chr1	15	27	a4	4	-	1" > exp
$BT multicov -bams two_blocks.bam -bed multicov.bed -split -f 0.01 > obs
check obs exp
rm obs exp

##################################################################
#  Test split alignment matches with -split and -f
##################################################################
echo -e "    multicov.t9...\c"
echo \
"chr1	15	20	a1	1	+	0
chr1	15	27	a2	2	+	0
chr1	15	20	a3	3	-	0
chr1	15	27	a4	4	-	0" > exp
$BT multicov -bams two_blocks.bam -bed multicov.bed -split -f 0.10 > obs
check obs exp
rm obs exp


##################################################################
#  Test when one of the iterator returns an empty set
##################################################################
echo -e "    multicov.t10...\c"
echo \
"chr1	0	250	4	0
chr1	500	1000	0	4"	>	exp
$BT multicov -bams test-multi.bam test-multi.2.bam -bed test-multi.bed > obs
check obs exp
rm obs exp



rm *.bam *.bai
[[ $FAILURES -eq 0 ]] || exit 1;
