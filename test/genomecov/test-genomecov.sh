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
###########################################################
#                       BAM files                         #
###########################################################
###########################################################
samtools view -Sb one_block.sam > one_block.bam 2>/dev/null
samtools view -Sb two_blocks.sam > two_blocks.bam 2>/dev/null
samtools view -Sb three_blocks.sam > three_blocks.bam 2>/dev/null
samtools view -Sb sam-w-del.sam > sam-w-del.bam 2>/dev/null



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


rm *.bam