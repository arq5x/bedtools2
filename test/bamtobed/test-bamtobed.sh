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

###########################################################
###########################################################
#                       BAM files                         #
###########################################################
###########################################################
$htsutil samtobam one_block.sam one_block.bam
$htsutil samtobam two_blocks.sam two_blocks.bam
$htsutil samtobam three_blocks.sam three_blocks.bam
#$htsutil samtobam sam-w-del.sam sam-w-del.bam
$htsutil samtobam two_blocks_w_D.sam two_blocks_w_D.bam
$htsutil samtobam numeric_tag.sam numeric_tag.bam
gzip -dfc extra-long-header.sam | $htsutil samtobam - extra-long-header.bam


##################################################################
#  Test one block without -split
##################################################################
echo -e "    bamtobed.t1...\c"
echo \
"chr1	0	30	one_blocks	40	-" > exp
$BT bamtobed -i one_block.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test one block with -split
##################################################################
echo -e "    bamtobed.t2...\c"
echo \
"chr1	0	30	one_blocks	40	-" > exp
$BT bamtobed -i one_block.bam -split > obs
check obs exp
rm obs exp


##################################################################
#  Test two blocks without -split
##################################################################
echo -e "    bamtobed.t3...\c"
echo \
"chr1	0	40	two_blocks	40	-" > exp
$BT bamtobed -i two_blocks.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test two blocks with -split
##################################################################
echo -e "    bamtobed.t4...\c"
echo \
"chr1	0	15	two_blocks	40	-
chr1	25	40	two_blocks	40	-" > exp
$BT bamtobed -i two_blocks.bam -split > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks without -split
##################################################################
echo -e "    bamtobed.t5...\c"
echo \
"chr1	0	50	three_blocks	40	-" > exp
$BT bamtobed -i three_blocks.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test three blocks with -split
##################################################################
echo -e "    bamtobed.t6...\c"
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
echo -e "    bamtobed.t7...\c"
echo \
"chr1	0	50	three_blocks	40	-	0	50	255,0,0	3	10,10,10	0,20,40" > exp
$BT bamtobed -i three_blocks.bam -bed12 > obs
check obs exp
rm obs exp


##################################################################
#  Ensure that both ways of getting blocks from a spliced alignment
# are indenticsl
##################################################################
echo -e "    bamtobed.t8...\c"
$BT bamtobed -i three_blocks.bam -split > split
$BT bamtobed -i three_blocks.bam -bed12 |  $BT bed12tobed6 > bed12
check split bed12
rm split bed12

##################################################################
# Test an alignment with a D operator and N operator -split option
##################################################################
echo -e "    bamtobed.t9...\c"
echo \
"chr1	0	15	two_blocks_1_1/2	40	+
chr1	25	40	two_blocks_1_1/2	40	+
chr1	99	129	two_blocks_1_2/1	40	+
chr1	0	15	two_blocks_2_1/2	40	+
chr1	25	42	two_blocks_2_1/2	40	+
chr1	99	129	two_blocks_2_2/1	40	+" > exp
$BT bamtobed -i two_blocks_w_D.bam -split > obs
check exp obs
rm exp obs


##################################################################
# Test an alignment with a D operator and N operator -splitD option
##################################################################
echo -e "    bamtobed.t10...\c"
echo \
"chr1	0	15	two_blocks_1_1/2	40	+
chr1	25	40	two_blocks_1_1/2	40	+
chr1	99	129	two_blocks_1_2/1	40	+
chr1	0	15	two_blocks_2_1/2	40	+
chr1	25	35	two_blocks_2_1/2	40	+
chr1	37	42	two_blocks_2_1/2	40	+
chr1	99	129	two_blocks_2_2/1	40	+" > exp
$BT bamtobed -i two_blocks_w_D.bam -splitD > obs
check exp obs
rm exp obs


##################################################################
# Test an alignment with a D operator and N operator -bed12 option
##################################################################
echo -e "    bamtobed.t9...\c"
echo \
"chr1	0	40	two_blocks_1_1/2	40	+	0	40	255,0,0	2	15,15	0,25
chr1	99	129	two_blocks_1_2/1	40	+	99	129	255,0,0	1	30	0
chr1	0	42	two_blocks_2_1/2	40	+	0	42	255,0,0	2	15,17	0,25
chr1	99	129	two_blocks_2_2/1	40	+	99	129	255,0,0	1	30	0" > exp
$BT bamtobed -i two_blocks_w_D.bam -bed12 > obs
check exp obs
rm exp obs


##################################################################
# Test an alignment with a D operator and N operator -bed12 option
##################################################################
echo -e "    bamtobed.t11...\c"
echo \
"chr1	0	40	two_blocks_1_1/2	40	+	0	40	255,0,0	2	15,15	0,25
chr1	99	129	two_blocks_1_2/1	40	+	99	129	255,0,0	1	30	0
chr1	0	42	two_blocks_2_1/2	40	+	0	42	255,0,0	3	15,10,5	0,25,37
chr1	99	129	two_blocks_2_2/1	40	+	99	129	255,0,0	1	30	0" > exp
$BT bamtobed -i two_blocks_w_D.bam -bed12 -splitD > obs

check exp obs
rm exp obs


##################################################################
# Test BAM input with numeric Tags
##################################################################
echo -e "    bamtobed.t12..\c"
echo \
"1	9998	10056	HISEQ1:18:H8VC6ADXX:1:1201:3360:80789/2	2	+
1	9998	10064	HISEQ1:18:H8VC6ADXX:1:1208:16920:47717/2	2	+"	>	exp
$BT bamtobed -tag NM -i numeric_tag.bam > obs

check exp obs
rm exp obs


##################################################################
# Test BAM input with long header
##################################################################
echo -e "    bamtobed.t13..\c"
echo -n ""	>	exp
$BT bamtobed -tag NM -i extra-long-header.bam > obs
check exp obs
rm exp obs


rm *.bam


[[ $FAILURES -eq 0 ]] || exit 1;
