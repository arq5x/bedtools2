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
$htsutil samtobam sam-w-del.sam sam-w-del.bam


##################################################################
#  Test three blocks without -split
##################################################################
echo -e "    coverage.t1...\c"
echo \
"chr1	0	50	three_blocks	40	-	0	50	0,0,0	3	10,10,10,	0,20,40,	1	40	50	0.8000000" > exp
$BT coverage -abam three_blocks.bam -b three_blocks_nomatch.bed > obs
check obs exp
echo -e "    coverage.t1b...\c"
$BT coverage -abam three_blocks.bam -b three_blocks_nomatch.bed -sorted > obs
check obs exp
rm obs exp


##################################################################
#  Test default
##################################################################
echo -e "    coverage.t2...\c"
echo \
"chr1	20	70	6	25	+	2	50	50	1.0000000
chr1	50	100	1	25	-	5	50	50	1.0000000
chr1	200	250	3	25	+	4	38	50	0.7600000
chr2	80	130	5	25	-	6	50	50	1.0000000
chr2	150	200	4	25	+	7	50	50	1.0000000
chr2	180	230	2	25	-	6	50	50	1.0000000" > exp
$BT coverage -a a.bed -b b.bed > obs
check exp obs
echo -e "    coverage.t2b...\c"
$BT coverage -a a.bed -b b.bed -sorted > obs
check exp obs
rm exp obs


##################################################################
#  Test counts
##################################################################
echo -e "    coverage.t3...\c"
echo \
"chr1	20	70	6	25	+	2
chr1	50	100	1	25	-	5
chr1	200	250	3	25	+	4
chr2	80	130	5	25	-	6
chr2	150	200	4	25	+	7
chr2	180	230	2	25	-	6" > exp
$BT coverage -a a.bed -b b.bed -counts > obs
check exp obs
echo -e "    coverage.t3b...\c"
$BT coverage -a a.bed -b b.bed -counts -sorted > obs
check exp obs
rm exp obs



##################################################################
#  Test -hist
##################################################################
echo -e "    coverage.t4...\c"
echo \
"chr1	20	70	6	25	+	2	50	50	1.0000000
chr1	50	100	1	25	-	2	40	50	0.8000000
chr1	50	100	1	25	-	3	10	50	0.2000000
chr1	200	250	3	25	+	0	12	50	0.2400000
chr1	200	250	3	25	+	1	20	50	0.4000000
chr1	200	250	3	25	+	2	11	50	0.2200000
chr1	200	250	3	25	+	3	4	50	0.0800000
chr1	200	250	3	25	+	4	3	50	0.0600000
chr2	80	130	5	25	-	3	46	50	0.9200000
chr2	80	130	5	25	-	4	4	50	0.0800000
chr2	150	200	4	25	+	5	22	50	0.4400000
chr2	150	200	4	25	+	6	28	50	0.5600000
chr2	180	230	2	25	-	1	16	50	0.3200000
chr2	180	230	2	25	-	2	2	50	0.0400000
chr2	180	230	2	25	-	3	6	50	0.1200000
chr2	180	230	2	25	-	4	4	50	0.0800000
chr2	180	230	2	25	-	5	13	50	0.2600000
chr2	180	230	2	25	-	6	9	50	0.1800000
all	0	12	300	0.0400000
all	1	36	300	0.1200000
all	2	103	300	0.3433333
all	3	66	300	0.2200000
all	4	11	300	0.0366667
all	5	35	300	0.1166667
all	6	37	300	0.1233333" > exp
$BT coverage -a a.bed -b b.bed -hist > obs
check exp obs
echo -e "    coverage.t4b...\c"
$BT coverage -a a.bed -b b.bed -hist -sorted > obs
check exp obs
rm exp obs

##################################################################
#  Test -d
##################################################################
echo -e "    coverage.t5...\c"
echo \
"chr1	20	70	6	25	+	1	2
chr1	20	70	6	25	+	2	2
chr1	20	70	6	25	+	3	2
chr1	20	70	6	25	+	4	2
chr1	20	70	6	25	+	5	2
chr1	20	70	6	25	+	6	2
chr1	20	70	6	25	+	7	2
chr1	20	70	6	25	+	8	2
chr1	20	70	6	25	+	9	2
chr1	20	70	6	25	+	10	2
chr1	20	70	6	25	+	11	2
chr1	20	70	6	25	+	12	2
chr1	20	70	6	25	+	13	2
chr1	20	70	6	25	+	14	2
chr1	20	70	6	25	+	15	2
chr1	20	70	6	25	+	16	2
chr1	20	70	6	25	+	17	2
chr1	20	70	6	25	+	18	2
chr1	20	70	6	25	+	19	2
chr1	20	70	6	25	+	20	2
chr1	20	70	6	25	+	21	2
chr1	20	70	6	25	+	22	2
chr1	20	70	6	25	+	23	2
chr1	20	70	6	25	+	24	2
chr1	20	70	6	25	+	25	2
chr1	20	70	6	25	+	26	2
chr1	20	70	6	25	+	27	2
chr1	20	70	6	25	+	28	2
chr1	20	70	6	25	+	29	2
chr1	20	70	6	25	+	30	2
chr1	20	70	6	25	+	31	2
chr1	20	70	6	25	+	32	2
chr1	20	70	6	25	+	33	2
chr1	20	70	6	25	+	34	2
chr1	20	70	6	25	+	35	2
chr1	20	70	6	25	+	36	2
chr1	20	70	6	25	+	37	2
chr1	20	70	6	25	+	38	2
chr1	20	70	6	25	+	39	2
chr1	20	70	6	25	+	40	2
chr1	20	70	6	25	+	41	2
chr1	20	70	6	25	+	42	2
chr1	20	70	6	25	+	43	2
chr1	20	70	6	25	+	44	2
chr1	20	70	6	25	+	45	2
chr1	20	70	6	25	+	46	2
chr1	20	70	6	25	+	47	2
chr1	20	70	6	25	+	48	2
chr1	20	70	6	25	+	49	2
chr1	20	70	6	25	+	50	2
chr1	50	100	1	25	-	1	2
chr1	50	100	1	25	-	2	2
chr1	50	100	1	25	-	3	2
chr1	50	100	1	25	-	4	2
chr1	50	100	1	25	-	5	2
chr1	50	100	1	25	-	6	2
chr1	50	100	1	25	-	7	2
chr1	50	100	1	25	-	8	2
chr1	50	100	1	25	-	9	2
chr1	50	100	1	25	-	10	2
chr1	50	100	1	25	-	11	2
chr1	50	100	1	25	-	12	2
chr1	50	100	1	25	-	13	2
chr1	50	100	1	25	-	14	2
chr1	50	100	1	25	-	15	2
chr1	50	100	1	25	-	16	2
chr1	50	100	1	25	-	17	2
chr1	50	100	1	25	-	18	2
chr1	50	100	1	25	-	19	2
chr1	50	100	1	25	-	20	2
chr1	50	100	1	25	-	21	2
chr1	50	100	1	25	-	22	2
chr1	50	100	1	25	-	23	2
chr1	50	100	1	25	-	24	2
chr1	50	100	1	25	-	25	3
chr1	50	100	1	25	-	26	2
chr1	50	100	1	25	-	27	2
chr1	50	100	1	25	-	28	3
chr1	50	100	1	25	-	29	2
chr1	50	100	1	25	-	30	2
chr1	50	100	1	25	-	31	2
chr1	50	100	1	25	-	32	2
chr1	50	100	1	25	-	33	2
chr1	50	100	1	25	-	34	2
chr1	50	100	1	25	-	35	2
chr1	50	100	1	25	-	36	2
chr1	50	100	1	25	-	37	2
chr1	50	100	1	25	-	38	2
chr1	50	100	1	25	-	39	2
chr1	50	100	1	25	-	40	2
chr1	50	100	1	25	-	41	2
chr1	50	100	1	25	-	42	2
chr1	50	100	1	25	-	43	3
chr1	50	100	1	25	-	44	3
chr1	50	100	1	25	-	45	3
chr1	50	100	1	25	-	46	3
chr1	50	100	1	25	-	47	3
chr1	50	100	1	25	-	48	3
chr1	50	100	1	25	-	49	3
chr1	50	100	1	25	-	50	3
chr1	200	250	3	25	+	1	4
chr1	200	250	3	25	+	2	4
chr1	200	250	3	25	+	3	4
chr1	200	250	3	25	+	4	3
chr1	200	250	3	25	+	5	3
chr1	200	250	3	25	+	6	3
chr1	200	250	3	25	+	7	3
chr1	200	250	3	25	+	8	2
chr1	200	250	3	25	+	9	2
chr1	200	250	3	25	+	10	2
chr1	200	250	3	25	+	11	2
chr1	200	250	3	25	+	12	2
chr1	200	250	3	25	+	13	2
chr1	200	250	3	25	+	14	2
chr1	200	250	3	25	+	15	2
chr1	200	250	3	25	+	16	2
chr1	200	250	3	25	+	17	2
chr1	200	250	3	25	+	18	2
chr1	200	250	3	25	+	19	1
chr1	200	250	3	25	+	20	1
chr1	200	250	3	25	+	21	1
chr1	200	250	3	25	+	22	1
chr1	200	250	3	25	+	23	1
chr1	200	250	3	25	+	24	1
chr1	200	250	3	25	+	25	1
chr1	200	250	3	25	+	26	1
chr1	200	250	3	25	+	27	1
chr1	200	250	3	25	+	28	1
chr1	200	250	3	25	+	29	1
chr1	200	250	3	25	+	30	1
chr1	200	250	3	25	+	31	1
chr1	200	250	3	25	+	32	1
chr1	200	250	3	25	+	33	1
chr1	200	250	3	25	+	34	1
chr1	200	250	3	25	+	35	1
chr1	200	250	3	25	+	36	1
chr1	200	250	3	25	+	37	1
chr1	200	250	3	25	+	38	1
chr1	200	250	3	25	+	39	0
chr1	200	250	3	25	+	40	0
chr1	200	250	3	25	+	41	0
chr1	200	250	3	25	+	42	0
chr1	200	250	3	25	+	43	0
chr1	200	250	3	25	+	44	0
chr1	200	250	3	25	+	45	0
chr1	200	250	3	25	+	46	0
chr1	200	250	3	25	+	47	0
chr1	200	250	3	25	+	48	0
chr1	200	250	3	25	+	49	0
chr1	200	250	3	25	+	50	0
chr2	80	130	5	25	-	1	4
chr2	80	130	5	25	-	2	3
chr2	80	130	5	25	-	3	3
chr2	80	130	5	25	-	4	3
chr2	80	130	5	25	-	5	3
chr2	80	130	5	25	-	6	3
chr2	80	130	5	25	-	7	3
chr2	80	130	5	25	-	8	3
chr2	80	130	5	25	-	9	3
chr2	80	130	5	25	-	10	3
chr2	80	130	5	25	-	11	3
chr2	80	130	5	25	-	12	3
chr2	80	130	5	25	-	13	3
chr2	80	130	5	25	-	14	3
chr2	80	130	5	25	-	15	3
chr2	80	130	5	25	-	16	3
chr2	80	130	5	25	-	17	3
chr2	80	130	5	25	-	18	3
chr2	80	130	5	25	-	19	3
chr2	80	130	5	25	-	20	3
chr2	80	130	5	25	-	21	3
chr2	80	130	5	25	-	22	3
chr2	80	130	5	25	-	23	3
chr2	80	130	5	25	-	24	3
chr2	80	130	5	25	-	25	3
chr2	80	130	5	25	-	26	3
chr2	80	130	5	25	-	27	3
chr2	80	130	5	25	-	28	3
chr2	80	130	5	25	-	29	3
chr2	80	130	5	25	-	30	3
chr2	80	130	5	25	-	31	3
chr2	80	130	5	25	-	32	3
chr2	80	130	5	25	-	33	3
chr2	80	130	5	25	-	34	3
chr2	80	130	5	25	-	35	3
chr2	80	130	5	25	-	36	3
chr2	80	130	5	25	-	37	3
chr2	80	130	5	25	-	38	3
chr2	80	130	5	25	-	39	3
chr2	80	130	5	25	-	40	3
chr2	80	130	5	25	-	41	3
chr2	80	130	5	25	-	42	3
chr2	80	130	5	25	-	43	3
chr2	80	130	5	25	-	44	3
chr2	80	130	5	25	-	45	3
chr2	80	130	5	25	-	46	3
chr2	80	130	5	25	-	47	3
chr2	80	130	5	25	-	48	4
chr2	80	130	5	25	-	49	4
chr2	80	130	5	25	-	50	4
chr2	150	200	4	25	+	1	6
chr2	150	200	4	25	+	2	6
chr2	150	200	4	25	+	3	5
chr2	150	200	4	25	+	4	5
chr2	150	200	4	25	+	5	5
chr2	150	200	4	25	+	6	5
chr2	150	200	4	25	+	7	5
chr2	150	200	4	25	+	8	5
chr2	150	200	4	25	+	9	5
chr2	150	200	4	25	+	10	5
chr2	150	200	4	25	+	11	5
chr2	150	200	4	25	+	12	5
chr2	150	200	4	25	+	13	5
chr2	150	200	4	25	+	14	6
chr2	150	200	4	25	+	15	6
chr2	150	200	4	25	+	16	6
chr2	150	200	4	25	+	17	6
chr2	150	200	4	25	+	18	6
chr2	150	200	4	25	+	19	6
chr2	150	200	4	25	+	20	6
chr2	150	200	4	25	+	21	6
chr2	150	200	4	25	+	22	6
chr2	150	200	4	25	+	23	6
chr2	150	200	4	25	+	24	6
chr2	150	200	4	25	+	25	6
chr2	150	200	4	25	+	26	6
chr2	150	200	4	25	+	27	6
chr2	150	200	4	25	+	28	6
chr2	150	200	4	25	+	29	6
chr2	150	200	4	25	+	30	6
chr2	150	200	4	25	+	31	6
chr2	150	200	4	25	+	32	6
chr2	150	200	4	25	+	33	6
chr2	150	200	4	25	+	34	6
chr2	150	200	4	25	+	35	6
chr2	150	200	4	25	+	36	6
chr2	150	200	4	25	+	37	6
chr2	150	200	4	25	+	38	6
chr2	150	200	4	25	+	39	6
chr2	150	200	4	25	+	40	5
chr2	150	200	4	25	+	41	5
chr2	150	200	4	25	+	42	5
chr2	150	200	4	25	+	43	5
chr2	150	200	4	25	+	44	5
chr2	150	200	4	25	+	45	5
chr2	150	200	4	25	+	46	5
chr2	150	200	4	25	+	47	5
chr2	150	200	4	25	+	48	5
chr2	150	200	4	25	+	49	5
chr2	150	200	4	25	+	50	5
chr2	180	230	2	25	-	1	6
chr2	180	230	2	25	-	2	6
chr2	180	230	2	25	-	3	6
chr2	180	230	2	25	-	4	6
chr2	180	230	2	25	-	5	6
chr2	180	230	2	25	-	6	6
chr2	180	230	2	25	-	7	6
chr2	180	230	2	25	-	8	6
chr2	180	230	2	25	-	9	6
chr2	180	230	2	25	-	10	5
chr2	180	230	2	25	-	11	5
chr2	180	230	2	25	-	12	5
chr2	180	230	2	25	-	13	5
chr2	180	230	2	25	-	14	5
chr2	180	230	2	25	-	15	5
chr2	180	230	2	25	-	16	5
chr2	180	230	2	25	-	17	5
chr2	180	230	2	25	-	18	5
chr2	180	230	2	25	-	19	5
chr2	180	230	2	25	-	20	5
chr2	180	230	2	25	-	21	5
chr2	180	230	2	25	-	22	5
chr2	180	230	2	25	-	23	4
chr2	180	230	2	25	-	24	4
chr2	180	230	2	25	-	25	4
chr2	180	230	2	25	-	26	4
chr2	180	230	2	25	-	27	3
chr2	180	230	2	25	-	28	3
chr2	180	230	2	25	-	29	3
chr2	180	230	2	25	-	30	3
chr2	180	230	2	25	-	31	3
chr2	180	230	2	25	-	32	3
chr2	180	230	2	25	-	33	2
chr2	180	230	2	25	-	34	2
chr2	180	230	2	25	-	35	1
chr2	180	230	2	25	-	36	1
chr2	180	230	2	25	-	37	1
chr2	180	230	2	25	-	38	1
chr2	180	230	2	25	-	39	1
chr2	180	230	2	25	-	40	1
chr2	180	230	2	25	-	41	1
chr2	180	230	2	25	-	42	1
chr2	180	230	2	25	-	43	1
chr2	180	230	2	25	-	44	1
chr2	180	230	2	25	-	45	1
chr2	180	230	2	25	-	46	1
chr2	180	230	2	25	-	47	1
chr2	180	230	2	25	-	48	1
chr2	180	230	2	25	-	49	1
chr2	180	230	2	25	-	50	1" > exp
$BT coverage -a a.bed -b b.bed -d > obs
check exp obs
echo -e "    coverage.t5b...\c"
$BT coverage -a a.bed -b b.bed -d -sorted > obs
check exp obs
rm exp obs


##################################################################
#  Test mean
##################################################################
echo -e "    coverage.t6...\c"
echo \
"chr1	20	70	6	25	+	2.0000000
chr1	50	100	1	25	-	2.2000000
chr1	200	250	3	25	+	1.3200001
chr2	80	130	5	25	-	3.0799999
chr2	150	200	4	25	+	5.5599999
chr2	180	230	2	25	-	3.4600000" > exp
$BT coverage -a a.bed -b b.bed -mean > obs
check exp obs
echo -e "    coverage.t6b...\c"
$BT coverage -a a.bed -b b.bed -mean -sorted > obs
check exp obs
rm exp obs


##################################################################
#  Test -s
##################################################################
echo -e "    coverage.t7...\c"
echo \
"chr1	20	70	6	25	+	2	50	50	1.0000000
chr1	50	100	1	25	-	1	23	50	0.4600000
chr1	200	250	3	25	+	0	0	50	0.0000000
chr2	80	130	5	25	-	4	50	50	1.0000000
chr2	150	200	4	25	+	3	50	50	1.0000000
chr2	180	230	2	25	-	4	34	50	0.6800000" > exp
$BT coverage -a a.bed -b b.bed -s > obs
check exp obs
echo -e "    coverage.t7b...\c"
$BT coverage -a a.bed -b b.bed -s -sorted > obs
check exp obs
rm exp obs


##################################################################
#  Test -S
##################################################################
echo -e "    coverage.t8...\c"
echo \
"chr1	20	70	6	25	+	0	0	50	0.0000000
chr1	50	100	1	25	-	4	50	50	1.0000000
chr1	200	250	3	25	+	4	38	50	0.7600000
chr2	80	130	5	25	-	2	50	50	1.0000000
chr2	150	200	4	25	+	4	50	50	1.0000000
chr2	180	230	2	25	-	2	50	50	1.0000000" > exp
$BT coverage -a a.bed -b b.bed -S > obs
check exp obs
echo -e "    coverage.t8b...\c"
$BT coverage -a a.bed -b b.bed -S -sorted > obs
check exp obs
rm exp obs

##################################################################
#  Test -S
##################################################################
echo -e "    coverage.t9...\c"
echo \
"chr1	20	70	6	25	+	0	0	50	0.0000000
chr1	50	100	1	25	-	4	50	50	1.0000000
chr1	200	250	3	25	+	4	38	50	0.7600000
chr2	80	130	5	25	-	2	50	50	1.0000000
chr2	150	200	4	25	+	4	50	50	1.0000000
chr2	180	230	2	25	-	2	50	50	1.0000000" > exp
$BT coverage -a a.bed -b b.bed -S > obs
check exp obs
echo -e "    coverage.t9b...\c"
$BT coverage -a a.bed -b b.bed -S -sorted > obs
check exp obs
rm exp obs


##################################################################
#  Test -split
##################################################################
echo -e "    coverage.t10...\c"
echo \
"chr1	0	50	3	30	50	0.6000000
chr1	12	20	0	0	8	0.0000000" > exp
$BT coverage -a c.bed -b three_blocks_match.bam -split > obs
check exp obs
echo -e "    coverage.t10b...\c"
$BT coverage -a c.bed -b three_blocks_match.bam -split -sorted > obs
check exp obs
rm exp obs

##################################################################
#  Test w/o -split
##################################################################
echo -e "    coverage.t11...\c"
echo \
"chr1	0	50	1	50	50	1.0000000
chr1	12	20	1	8	8	1.0000000" > exp
$BT coverage -a c.bed -b three_blocks_match.bam  > obs
check exp obs
echo -e "    coverage.t11b...\c"
$BT coverage -a c.bed -b three_blocks_match.bam -sorted > obs
check exp obs
rm exp obs

##################################################################
#  Test -split and -d
##################################################################
echo -e "    coverage.t12...\c"
echo \
"chr1	0	50	1	1
chr1	0	50	2	1
chr1	0	50	3	1
chr1	0	50	4	1
chr1	0	50	5	1
chr1	0	50	6	1
chr1	0	50	7	1
chr1	0	50	8	1
chr1	0	50	9	1
chr1	0	50	10	1
chr1	0	50	11	0
chr1	0	50	12	0
chr1	0	50	13	0
chr1	0	50	14	0
chr1	0	50	15	0
chr1	0	50	16	0
chr1	0	50	17	0
chr1	0	50	18	0
chr1	0	50	19	0
chr1	0	50	20	0
chr1	0	50	21	1
chr1	0	50	22	1
chr1	0	50	23	1
chr1	0	50	24	1
chr1	0	50	25	1
chr1	0	50	26	1
chr1	0	50	27	1
chr1	0	50	28	1
chr1	0	50	29	1
chr1	0	50	30	1
chr1	0	50	31	0
chr1	0	50	32	0
chr1	0	50	33	0
chr1	0	50	34	0
chr1	0	50	35	0
chr1	0	50	36	0
chr1	0	50	37	0
chr1	0	50	38	0
chr1	0	50	39	0
chr1	0	50	40	0
chr1	0	50	41	1
chr1	0	50	42	1
chr1	0	50	43	1
chr1	0	50	44	1
chr1	0	50	45	1
chr1	0	50	46	1
chr1	0	50	47	1
chr1	0	50	48	1
chr1	0	50	49	1
chr1	0	50	50	1
chr1	12	20	1	0
chr1	12	20	2	0
chr1	12	20	3	0
chr1	12	20	4	0
chr1	12	20	5	0
chr1	12	20	6	0
chr1	12	20	7	0
chr1	12	20	8	0" > exp
$BT coverage -a c.bed -b three_blocks_match.bam -split -d > obs
check exp obs
echo -e "    coverage.t12b...\c"
$BT coverage -a c.bed -b three_blocks_match.bam -split -d -sorted > obs
check exp obs
rm exp obs

##################################################################
#  Test -split and -hist
##################################################################
echo -e "    coverage.t13...\c"
echo \
"chr1	0	50	0	20	50	0.4000000
chr1	0	50	1	30	50	0.6000000
chr1	12	20	0	8	8	1.0000000
all	0	28	58	0.4827586
all	1	30	58	0.5172414" > exp
$BT coverage -a c.bed -b three_blocks_match.bam -split -hist > obs
check exp obs
echo -e "    coverage.t13b...\c"
$BT coverage -a c.bed -b three_blocks_match.bam -split -hist -sorted > obs
check exp obs
rm exp obs


##################################################################
#  Test that -counts, -hist are mutually exclusive options
##################################################################
echo -e "    coverage.t14...\c"
echo \
"***** ERROR: -counts, -d, -mean, and -hist are all mutually exclusive options. *****" > exp
$BT coverage -a a.bed -b b.bed -counts -hist 2>&1 > /dev/null | tail -1 | cat - > obs
check exp obs
echo -e "    coverage.t14b...\c"
$BT coverage -a a.bed -b b.bed -counts -sorted -hist 2>&1 > /dev/null | tail -1 | cat - > obs
check exp obs
rm exp obs

##################################################################
#  Test that -counts, -d are mutually exclusive options
##################################################################
echo -e "    coverage.t15...\c"
echo \
"***** ERROR: -counts, -d, -mean, and -hist are all mutually exclusive options. *****" > exp
$BT coverage -a a.bed -b b.bed -counts -d 2>&1 > /dev/null | tail -1 | cat - > obs
check exp obs
echo -e "    coverage.t15b...\c"
$BT coverage -a a.bed -b b.bed -sorted -counts -d 2>&1 > /dev/null | tail -1 | cat - > obs
check exp obs
rm exp obs

##################################################################
#  Test that -hist, -d are mutually exclusive options
##################################################################
echo -e "    coverage.t16...\c"
echo \
"***** ERROR: -counts, -d, -mean, and -hist are all mutually exclusive options. *****" > exp
$BT coverage -a a.bed -b b.bed -hist -d 2>&1 > /dev/null | tail -1 | cat - > obs
check exp obs
echo -e "    coverage.t16b...\c"
$BT coverage -a a.bed -b b.bed -sorted -hist -d 2>&1 > /dev/null | tail -1 | cat - > obs
check exp obs
rm exp obs


##################################################################
#  Test that -mean, -d are mutually exclusive options
##################################################################
echo -e "    coverage.t17...\c"
echo \
"***** ERROR: -counts, -d, -mean, and -hist are all mutually exclusive options. *****" > exp
$BT coverage -a a.bed -b b.bed -mean -d 2>&1 > /dev/null | tail -1 | cat - > obs
check exp obs
echo -e "    coverage.t17b...\c"
$BT coverage -a a.bed -b b.bed -sorted -mean -d 2>&1 > /dev/null | tail -1 | cat - > obs
check exp obs
rm exp obs

##################################################################
#  Test the last record in file with no overlaps is reported
##################################################################
echo -e "    coverage.t18...\c"
echo \
"chr1	0	10	1	0
chr1	0	10	2	0
chr1	0	10	3	0
chr1	0	10	4	1
chr1	0	10	5	1
chr1	0	10	6	1
chr1	0	10	7	1
chr1	0	10	8	1
chr1	0	10	9	1
chr1	0	10	10	1
chr1	15	20	1	0
chr1	15	20	2	0
chr1	15	20	3	0
chr1	15	20	4	0
chr1	15	20	5	0
chr1	21	25	1	0
chr1	21	25	2	0
chr1	21	25	3	0
chr1	21	25	4	0" > exp
$BT coverage -a x.bed -b y.bed -d > obs
check exp obs
echo -e "    coverage.t18b...\c"
$BT coverage -a x.bed -b y.bed -d -sorted > obs
check exp obs
rm exp obs

##################################################################
#  Test the last record in file with no overlaps is reported
##################################################################
echo -e "    coverage.t19...\c"
echo \
"chr1	0	10	0	3	10	0.3000000
chr1	0	10	1	7	10	0.7000000
chr1	15	20	0	5	5	1.0000000
chr1	21	25	0	4	4	1.0000000
all	0	12	19	0.6315789
all	1	7	19	0.3684210" > exp
$BT coverage -a x.bed -b y.bed -hist > obs
check exp obs
echo -e "    coverage.t19b...\c"
$BT coverage -a x.bed -b y.bed -hist -sorted > obs
check exp obs
rm exp obs

################################################################
# Test that simple chr	0	100 works
################################################################
echo -e "    coverage.t20...\c"
echo \
"chr	0	100	1	100	100	1.0000000" > exp
$BT coverage -a chr_0-100.bed -b chr_0-100.bed > obs
check exp obs
rm exp obs

rm one_block.bam two_blocks.bam three_blocks.bam sam-w-del.bam

[[ $FAILURES -eq 0 ]] || exit 1;
