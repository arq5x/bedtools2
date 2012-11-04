BT=${BT-../../bin/bedtools}

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

samtools view -Sb one_block.sam > one_block.bam 2>/dev/null
samtools view -Sb two_blocks.sam > two_blocks.bam 2>/dev/null
samtools index one_block.bam
samtools index two_blocks.bam

##################################################################
#  Test one block matches all BEDs
##################################################################
echo "    multicov.t1...\c"
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
echo "    multicov.t2...\c"
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
echo "    multicov.t3...\c"
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
echo "    multicov.t4...\c"
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
echo "    multicov.t5...\c"
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
echo "    multicov.t6...\c"
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
echo "    multicov.t7...\c"
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
echo "    multicov.t8...\c"
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
echo "    multicov.t9...\c"
echo \
"chr1	15	20	a1	1	+	0
chr1	15	27	a2	2	+	0
chr1	15	20	a3	3	-	0
chr1	15	27	a4	4	-	0" > exp
$BT multicov -bams two_blocks.bam -bed multicov.bed -split -f 0.10 > obs
check obs exp
rm obs exp



rm *.bam