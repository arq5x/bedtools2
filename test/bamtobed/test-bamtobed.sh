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
#  Test one block without -split
##################################################################
echo "    bamtobed.t1...\c"
echo \
"chr1	0	30	one_blocks	40	-" > exp
$BT bamtobed -i one_block.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test one block with -split
##################################################################
echo "    bamtobed.t2...\c"
echo \
"chr1	0	30	one_blocks	40	-" > exp
$BT bamtobed -i one_block.bam -split > obs
check obs exp
rm obs exp


##################################################################
#  Test two blocks without -split
##################################################################
echo "    bamtobed.t3...\c"
echo \
"chr1	0	40	two_blocks	40	-" > exp
$BT bamtobed -i two_blocks.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test two blocks with -split
##################################################################
echo "    bamtobed.t4...\c"
echo \
"chr1	0	15	two_blocks	40	-
chr1	25	40	two_blocks	40	-" > exp
$BT bamtobed -i two_blocks.bam -split > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks without -split
##################################################################
echo "    bamtobed.t5...\c"
echo \
"chr1	0	50	three_blocks	40	-" > exp
$BT bamtobed -i three_blocks.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks with -split
##################################################################
echo "    bamtobed.t6...\c"
echo \
"chr1	0	10	three_blocks	40	-
chr1	20	30	three_blocks	40	-
chr1	40	50	three_blocks	40	-" > exp
$BT bamtobed -i three_blocks.bam -split > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks with -bed12
##################################################################
echo "    bamtobed.t7...\c"
echo \
"chr1	0	50	three_blocks	40	-	0	50	255,0,0	3	10,10,10	0,20,40" > exp
$BT bamtobed -i three_blocks.bam -bed12 > obs
check obs exp
rm obs exp


##################################################################
#  Ensure that both ways of getting blocks from a spliced alignment
# are indenticsl
##################################################################
echo "    bamtobed.t8...\c"
$BT bamtobed -i three_blocks.bam -split > split
$BT bamtobed -i three_blocks.bam -bed12 |  $BT bed12tobed6 > bed12
check split bed12
rm split bed12

rm *.bam
