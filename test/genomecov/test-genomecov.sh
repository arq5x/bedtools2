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
###########################################################
#                       BAM files                         #
###########################################################
###########################################################
samtools view -Sb one_block.sam > one_block.bam 2>/dev/null
samtools view -Sb two_blocks.sam > two_blocks.bam 2>/dev/null
samtools view -Sb three_blocks.sam > three_blocks.bam 2>/dev/null
samtools view -Sb sam-w-del.sam > sam-w-del.bam 2>/dev/null
samtools view -Sb pair-chip.sam > pair-chip.bam 2>/dev/null
samtools view -Sb chip.sam > chip.bam 2>/dev/null



##################################################################
#  Test three blocks without -split
##################################################################
echo "    genomecov.t1...\c"
echo \
"chr1	0	50	1" > exp
$BT genomecov -ibam three_blocks.bam -bg > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks with -split
##################################################################
echo "    genomecov.t2...\c"
echo \
"chr1	0	10	1
chr1	20	30	1
chr1	40	50	1" > exp
$BT genomecov -ibam three_blocks.bam -bg -split > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks with -split and -bga
##################################################################
echo "    genomecov.t3...\c"
echo \
"chr1	0	10	1
chr1	10	20	0
chr1	20	30	1
chr1	30	40	0
chr1	40	50	1
chr1	50	1000	0" > exp
$BT genomecov -ibam three_blocks.bam -bga -split > obs
check obs exp
rm obs exp


##################################################################
#  Test blocked BAM from multiple files w/ -bga and w/o -split
##################################################################
echo "    genomecov.t4...\c"
echo \
"chr1	0	30	3
chr1	30	40	2
chr1	40	50	1
chr1	50	1000	0" > exp
samtools merge -f /dev/stdout *block*.bam | $BT genomecov -ibam - -bga  > obs
check obs exp
rm obs exp


##################################################################
#  Test blocked BAM from multiple files w/ -bga and w -split
##################################################################
echo "    genomecov.t5...\c"
echo \
"chr1	0	10	3
chr1	10	15	2
chr1	15	20	1
chr1	20	25	2
chr1	25	30	3
chr1	30	50	1
chr1	50	1000	0" > exp
samtools merge -f /dev/stdout *block*.bam | $BT genomecov -ibam - -bga -split > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks with -split and -dz
##################################################################
echo "    genomecov.t6...\c"
echo \
"chr1	0	1
chr1	1	1
chr1	2	1
chr1	3	1
chr1	4	1
chr1	5	1
chr1	6	1
chr1	7	1
chr1	8	1
chr1	9	1
chr1	20	1
chr1	21	1
chr1	22	1
chr1	23	1
chr1	24	1
chr1	25	1
chr1	26	1
chr1	27	1
chr1	28	1
chr1	29	1
chr1	40	1
chr1	41	1
chr1	42	1
chr1	43	1
chr1	44	1
chr1	45	1
chr1	46	1
chr1	47	1
chr1	48	1
chr1	49	1" > exp
$BT genomecov -ibam three_blocks.bam -dz -split > obs
check obs exp
rm obs exp


##################################################################
#  Test SAM with 1bp D operator
##################################################################
echo "    genomecov.t7...\c"
echo \
"chr1	0	10	1
chr1	11	21	1" > exp
$BT genomecov -ibam sam-w-del.bam -bg > obs
check obs exp
rm obs exp

##################################################################
#  Test bam with chroms that have no coverage
##################################################################
echo "    genomecov.t8...\c"
echo \
"1	0	93	100	0.93
1	1	4	100	0.04
1	2	3	100	0.03
2	0	100	100	1
3	0	100	100	1
genome	0	293	300	0.976667
genome	1	4	300	0.0133333
genome	2	3	300	0.01" > exp
$BT genomecov -ibam y.bam > obs
check obs exp
rm obs exp

##################################################################
#  Test bam with chroms that have no coverage
##################################################################
echo "    genomecov.t9...\c"
echo \
"1	15	17	1
1	17	20	2
1	20	22	1" > exp
$BT genomecov -ibam y.bam -bg > obs
check obs exp
rm obs exp

##################################################################
#  Test bam with chroms that have no coverage
##################################################################
echo "    genomecov.t10...\c"
echo \
"1	0	15	0
1	15	17	1
1	17	20	2
1	20	22	1
1	22	100	0
2	0	100	0
3	0	100	0" > exp
$BT genomecov -ibam y.bam -bga > obs
check obs exp
rm obs exp

##################################################################
#  Test bed with chroms that have no coverage
##################################################################
echo "    genomecov.t11...\c"
echo \
"1	0	93	100	0.93
1	1	4	100	0.04
1	2	3	100	0.03
2	0	100	100	1
3	0	100	100	1
genome	0	293	300	0.976667
genome	1	4	300	0.0133333
genome	2	3	300	0.01" > exp
$BT genomecov -i y.bed -g genome.txt > obs
check obs exp
rm obs exp

##################################################################
#  Test bed with chroms that have no coverage
##################################################################
echo "    genomecov.t12...\c"
echo \
"1	15	17	1
1	17	20	2
1	20	22	1" > exp
$BT genomecov -i y.bed -g genome.txt -bg > obs
check obs exp
rm obs exp

##################################################################
#  Test bed with chroms that have no coverage
##################################################################
echo "    genomecov.t13...\c"
echo \
"1	0	15	0
1	15	17	1
1	17	20	2
1	20	22	1
1	22	100	0
2	0	100	0
3	0	100	0" > exp
$BT genomecov -i y.bed -g genome.txt -bga > obs
check obs exp
rm obs exp

##################################################################
#  Test pair-end chip 
##################################################################
echo "    genomecov.t14...\c"
echo \
"chr1	0	203	1" > exp
$BT genomecov -ibam pair-chip.bam -bg -pc > obs
check obs exp
rm obs exp

##################################################################
#  Test chip fragmentSize
##################################################################
echo "    genomecov.t15...\c"
echo \
"chr1	1	101	1
chr1	200	300	1" > exp
$BT genomecov -ibam chip.bam -bg -fs 100 > obs
check obs exp
rm obs exp

