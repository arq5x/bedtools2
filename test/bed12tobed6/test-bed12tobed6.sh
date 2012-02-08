BT=../../bin/bedtools

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}


##################################################################
#  Test one block
##################################################################
echo "    bed12tobed6.t1...\c"
echo \
"chr1	0	50	one_blocks_match	0	+" > exp
$BT bed12tobed6 -i one_blocks.bed > obs
check obs exp
rm obs exp


##################################################################
#  Test two blocks
##################################################################
echo "    bed12tobed6.t2...\c"
echo \
"chr1	0	10	two_blocks_match	0	+
chr1	40	50	two_blocks_match	0	+" > exp
$BT bed12tobed6 -i two_blocks.bed > obs
check obs exp
rm obs exp

##################################################################
#  Test three blocks
##################################################################
echo "    bed12tobed6.t3...\c"
echo \
"chr1	0	10	three_blocks_match	0	+
chr1	20	30	three_blocks_match	0	+
chr1	40	50	three_blocks_match	0	+" > exp
$BT bed12tobed6 -i three_blocks.bed > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks and add block numbers
##################################################################
echo "    bed12tobed6.t4...\c"
echo \
"chr1	0	10	three_blocks_match	1	+
chr1	20	30	three_blocks_match	2	+
chr1	40	50	three_blocks_match	3	+" > exp
$BT bed12tobed6 -i three_blocks.bed -n > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks and add block numbers.  Test reverse strand
##################################################################
echo "    bed12tobed6.t5...\c"
echo \
"chr1	0	10	three_blocks_match	3	-
chr1	20	30	three_blocks_match	2	-
chr1	40	50	three_blocks_match	1	-" > exp
sed -e 's/\+/\-/' three_blocks.bed | $BT bed12tobed6 -n > obs
check obs exp
rm obs exp