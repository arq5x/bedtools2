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

###########################################################
###########################################################
#  coverage.new.*: tests added for the intersect rewrite.
#  coverage runs on intersect's readers, bin index and
#  overlap test (ContextCoverage subclasses ContextIntersect),
#  so each test pins a visible output of that machinery.
#  The fraction and mean columns are a float printed with
#  %0.7f; t002 and t031 pin that rounding.
###########################################################
###########################################################

###########################################################
#  coverage.new.t001: default columns for no, partial and
#  full-length hits
#  a1 (0-3) has B 0-1: 1 hit, 1 of 3 bases, 0.3333333.
#  a2 (10-13) has B 11-13: 1 hit, 2 of 3 bases, 0.6666667.
#  a3 (20-23) has no B: 0 0 3 0.0000000.
###########################################################
echo -e "    coverage.new.t001...\c"
printf 'chr1\t0\t3\ta1\nchr1\t10\t13\ta2\nchr1\t20\t23\ta3\n' > a.tmp.bed
printf 'chr1\t0\t1\nchr1\t11\t13\nchr1\t30\t40\n' > b.tmp.bed
printf 'chr1\t0\t3\ta1\t1\t1\t3\t0.3333333\nchr1\t10\t13\ta2\t1\t2\t3\t0.6666667\nchr1\t20\t23\ta3\t0\t0\t3\t0.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t002: the fraction is rounded as a float
#  7 of 19 bases covered: 7/19 = 0.368421052... in double,
#  but the float nearest 7/19 is 0.3684210479..., printed as
#  0.3684210. Pins the current output: float, not double.
###########################################################
echo -e "    coverage.new.t002...\c"
printf 'chr1\t0\t19\n' > a.tmp.bed
printf 'chr1\t5\t12\n' > b.tmp.bed
printf 'chr1\t0\t19\t1\t7\t19\t0.3684210\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t003: touching B are not counted, 1 bp is
#  B 200-300 and 50-100 only touch A 100-200 (half-open);
#  B 199-250 shares base 199: 1 hit, 1 base, 0.0100000.
###########################################################
echo -e "    coverage.new.t003...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t200\t300\nchr1\t50\t100\nchr1\t199\t250\n' > b.tmp.bed
printf 'chr1\t100\t200\t1\t1\t100\t0.0100000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t004: many overlapping B
#  Six B overlap A 100-200 (count 6). Their union inside A is
#  100-105, 110-140 and 160-175: 50 bases, 0.5000000.
###########################################################
echo -e "    coverage.new.t004...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t110\t130\nchr1\t120\t140\nchr1\t125\t128\nchr1\t160\t170\nchr1\t165\t175\nchr1\t90\t105\n' > b.tmp.bed
printf 'chr1\t100\t200\t6\t50\t100\t0.5000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t005: duplicated B are counted twice
#  Two copies of B 100-105 on A 100-110: count 2, 5 bases
#  covered (0.5000000); with -d bases 1-5 have depth 2.
###########################################################
echo -e "    coverage.new.t005...\c"
printf 'chr1\t100\t110\n' > a.tmp.bed
printf 'chr1\t100\t105\nchr1\t100\t105\n' > b.tmp.bed
printf 'chr1\t100\t110\t2\t5\t10\t0.5000000\n' > exp
for i in 1 2 3 4 5; do printf 'chr1\t100\t110\t%d\t2\n' $i >> exp; done
for i in 6 7 8 9 10; do printf 'chr1\t100\t110\t%d\t0\n' $i >> exp; done
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
$BT coverage -a a.tmp.bed -b b.tmp.bed -d >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t006: duplicated A records give two lines
###########################################################
echo -e "    coverage.new.t006...\c"
printf 'chr1\t100\t200\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\t1\t10\t100\t0.1000000\nchr1\t100\t200\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t007: unsorted A is reported in input order
#  q1 and q2 have no B; q3 has 10 of 100 bases covered; q4
#  (chr10:5-15) has 5 of 10 covered by chr10:0-10.
###########################################################
echo -e "    coverage.new.t007...\c"
printf 'chr2\t50\t60\tq1\nchr1\t300\t400\tq2\nchr1\t100\t200\tq3\nchr10\t5\t15\tq4\n' > a.tmp.bed
printf 'chr1\t150\t160\nchr10\t0\t10\n' > b.tmp.bed
printf 'chr2\t50\t60\tq1\t0\t0\t10\t0.0000000\nchr1\t300\t400\tq2\t0\t0\t100\t0.0000000\nchr1\t100\t200\tq3\t1\t10\t100\t0.1000000\nchr10\t5\t15\tq4\t1\t5\t10\t0.5000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t008: chromosomes present in only one file
#  chr1 is only in A (zeros); chr2 and chrY are only in B
#  (no lines); chrX has 10 of 100 bases covered.
###########################################################
echo -e "    coverage.new.t008...\c"
printf 'chr1\t100\t200\nchrX\t100\t200\n' > a.tmp.bed
printf 'chr2\t100\t200\nchrX\t150\t160\nchrY\t0\t1000\n' > b.tmp.bed
printf 'chr1\t100\t200\t0\t0\t100\t0.0000000\nchrX\t100\t200\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t009: records straddling bin boundaries
#  Bin edges at 16384, 131072, 1048576 (chr1), 2^29 (chr2's
#  B crosses it, so it is stored in a 2^32 bp bin) and 2^26
#  (chr3's B).
#  B 16383-16384 and 16384-16385 lie in fine bins 0 and 1,
#  on either side of 16384; e8 and e9 need both of them.
#  e1 16380-16390: B 16383-16385, 16383-16384, 16384-16385
#     cover 16383-16384 -> 3 2 10 0.2
#  e2 131070-131080: B 131071-131073 and B 131000-131100
#     -> 2 hits, all 10 bases
#  e3 1048570-1048580: B 1048575-1048577 -> 1 2 10 0.2
#  e4 16384-16390: B 16383-16385 and 16384-16385 cover base
#     16384 -> 2 1 6, 1/6 = 0.1666667
#  e5 131060-131072: both 131071-131073 and 131000-131100
#     -> 2 hits, all 12 bases
#  e6 chr2:300000000-300000100 inside chr2:0-600000000
#  e7 chr3:67108800-67108900: 67108860-67108870 -> 1 10 100
#  e8 16383-16390 (starts on the last base of fine bin 0):
#     all three 16383-16385 B -> 3 2 7, 2/7 = 0.2857143
#  e9 16380-16385 (ends on the first base of fine bin 1):
#     all three 16383-16385 B -> 3 2 5 0.4
###########################################################
echo -e "    coverage.new.t009...\c"
printf 'chr1\t16380\t16390\te1\nchr1\t131070\t131080\te2\nchr1\t1048570\t1048580\te3\nchr1\t16384\t16390\te4\nchr1\t131060\t131072\te5\nchr2\t300000000\t300000100\te6\nchr3\t67108800\t67108900\te7\nchr1\t16383\t16390\te8\nchr1\t16380\t16385\te9\n' > a.tmp.bed
printf 'chr1\t16383\t16385\nchr1\t131071\t131073\nchr1\t1048575\t1048577\nchr1\t131000\t131100\nchr2\t0\t600000000\nchr3\t67108860\t67108870\nchr1\t16383\t16384\nchr1\t16384\t16385\n' > b.tmp.bed
printf 'chr1\t16380\t16390\te1\t3\t2\t10\t0.2000000\nchr1\t131070\t131080\te2\t2\t10\t10\t1.0000000\nchr1\t1048570\t1048580\te3\t1\t2\t10\t0.2000000\nchr1\t16384\t16390\te4\t2\t1\t6\t0.1666667\nchr1\t131060\t131072\te5\t2\t12\t12\t1.0000000\nchr2\t300000000\t300000100\te6\t1\t100\t100\t1.0000000\nchr3\t67108800\t67108900\te7\t1\t10\t100\t0.1000000\nchr1\t16383\t16390\te8\t3\t2\t7\t0.2857143\nchr1\t16380\t16385\te9\t3\t2\t5\t0.4000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t010: 100 one-base B in and around one A
#  A chr1:16000-17000 straddles the 16384 bin edge. B has a
#  1 bp record every 10 bases from 15990 to 17010; the 100 at
#  16000, 16010, ..., 16990 fall in A: 100 hits, 100 bases.
###########################################################
echo -e "    coverage.new.t010...\c"
printf 'chr1\t16000\t17000\n' > a.tmp.bed
: > b.tmp.bed
i=15990
while [ $i -le 17010 ]; do
	printf 'chr1\t%d\t%d\n' $i $((i + 1)) >> b.tmp.bed
	i=$((i + 10))
done
printf 'chr1\t16000\t17000\t100\t100\t1000\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t011: coordinates above 2^31
#  A1 3000000000-3000000100 has 10 bases covered; A2
#  4000000000-4000000050 has 10 bases covered by a B that
#  starts before it. A3 4294967290-4294967300 lies inside B
#  4294967200-4294967400, which crosses 2^32 and so is stored
#  in the coarsest bin level: 1 hit, all 10 bases.
###########################################################
echo -e "    coverage.new.t011...\c"
printf 'chr1\t3000000000\t3000000100\nchr1\t4000000000\t4000000050\nchr1\t4294967290\t4294967300\n' > a.tmp.bed
printf 'chr1\t3000000050\t3000000060\nchr1\t3999999000\t4000000010\nchr1\t4294967200\t4294967400\n' > b.tmp.bed
printf 'chr1\t3000000000\t3000000100\t1\t10\t100\t0.1000000\nchr1\t4000000000\t4000000050\t1\t10\t50\t0.2000000\nchr1\t4294967290\t4294967300\t1\t10\t10\t1.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t012: A longer than 1024 bp, then shorter A
#  L1 (3000 bp) has bases 2990-2999 covered by B 2990-5005:
#  10/3000 = 0.0033333. S1 (5000-5010) has 5 bases from the
#  same B. L2 (2000 bp) has none. No depth carries over
#  between records.
###########################################################
echo -e "    coverage.new.t012...\c"
printf 'chr1\t0\t3000\tL1\nchr1\t5000\t5010\tS1\nchr1\t10000\t12000\tL2\n' > a.tmp.bed
printf 'chr1\t2990\t5005\n' > b.tmp.bed
printf 'chr1\t0\t3000\tL1\t1\t10\t3000\t0.0033333\nchr1\t5000\t5010\tS1\t1\t5\t10\t0.5000000\nchr1\t10000\t12000\tL2\t0\t0\t2000\t0.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t013: BED4 and BED5 A keep their columns
###########################################################
echo -e "    coverage.new.t013...\c"
printf 'chr1\t100\t200\tn1\n' > a4.tmp.bed
printf 'chr1\t100\t200\tn1\t9\n' > a5.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\tn1\t1\t10\t100\t0.1000000\nchr1\t100\t200\tn1\t9\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a4.tmp.bed -b b.tmp.bed > obs
$BT coverage -a a5.tmp.bed -b b.tmp.bed >> obs
check obs exp
rm obs exp a4.tmp.bed a5.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t014: BED6 plus extra columns keeps all of them
###########################################################
echo -e "    coverage.new.t014...\c"
printf 'chr1\t100\t200\tn1\t9\t+\textra1\textra2\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\tn1\t9\t+\textra1\textra2\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t015: bedGraph A keeps its value column
###########################################################
echo -e "    coverage.new.t015...\c"
printf 'chr1\t100\t200\t0.75\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\t0.75\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t016: BED12 A with and without -split
#  A spans 100-200 with blocks 100-120 and 180-200. B 130-150
#  is in the intron; B 190-195 is in block 2.
#  Without -split: 2 hits, 25 bases of the 100 bp span.
#  -counts -split: only 190-195 overlaps a block -> 1.
#  (The -split denominator is in bugs.md and is not asserted.)
###########################################################
echo -e "    coverage.new.t016...\c"
printf 'chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\n' > a.tmp.bed
printf 'chr1\t130\t150\nchr1\t190\t195\n' > b.tmp.bed
printf 'chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\t2\t25\t100\t0.2500000\nchr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\t1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
$BT coverage -a a.tmp.bed -b b.tmp.bed -counts -split >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t017: GFF A, default and -d
#  GFF 101-110 is 0-based 100-110 (10 bp). B 105-107 covers
#  0-based 105 and 106, which are positions 6 and 7 of A.
###########################################################
echo -e "    coverage.new.t017...\c"
printf 'chr1\tsrc\tgene\t101\t110\t.\t+\t.\tID=g1\n' > a.tmp.gff
printf 'chr1\t105\t107\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t101\t110\t.\t+\t.\tID=g1\t1\t2\t10\t0.2000000\n' > exp
for i in 1 2 3 4 5 6 7 8 9 10; do
	d=0
	if [ $i -eq 6 ] || [ $i -eq 7 ]; then d=1; fi
	printf 'chr1\tsrc\tgene\t101\t110\t.\t+\t.\tID=g1\t%d\t%d\n' $i $d >> exp
done
$BT coverage -a a.tmp.gff -b b.tmp.bed > obs
$BT coverage -a a.tmp.gff -b b.tmp.bed -d >> obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  coverage.new.t018: VCF A length is the REF length
#  rs1 (SNP, 0-based 124-125) and rs3 have no B. rs2 (REF
#  ACGT, 0-based 149-153) has B 148-152 over 149-151:
#  3 of 4 bases, 0.7500000.
###########################################################
echo -e "    coverage.new.t018...\c"
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t125\trs1\tA\tG\t50\tPASS\tDP=10\nchr1\t150\trs2\tACGT\tA\t50\tPASS\tDP=5\nchr1\t300\trs3\tC\tT\t50\tPASS\tDP=5\n' > a.tmp.vcf
printf 'chr1\t148\t152\n' > b.tmp.bed
printf 'chr1\t125\trs1\tA\tG\t50\tPASS\tDP=10\t0\t0\t1\t0.0000000\nchr1\t150\trs2\tACGT\tA\t50\tPASS\tDP=5\t1\t3\t4\t0.7500000\nchr1\t300\trs3\tC\tT\t50\tPASS\tDP=5\t0\t0\t1\t0.0000000\n' > exp
$BT coverage -a a.tmp.vcf -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  coverage.new.t019: BAM A, with and without -split
#  r1 is 50M at 100-150; r2 is 20M30N20M at 200-270 (blocks
#  200-220 and 250-270). Each read prints as BED12.
#  B 120-130 covers 10 bases of r1. B 225-245 lies wholly in
#  r2's N gap: without -split it covers 20 bases of r2's span
#  (20/70 = 0.2857143); -counts -split gives r1 1 and r2 0.
#  (The -split denominator is in bugs.md and is not asserted.)
###########################################################
echo -e "    coverage.new.t019...\c"
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*\nr2\t16\tchr1\t201\t30\t20M30N20M\t*\t0\t0\t*\t*\n' > a.tmp.sam
$htsutil samtobam a.tmp.sam a.tmp.bam
printf 'chr1\t120\t130\nchr1\t225\t245\n' > b.tmp.bed
printf 'chr1\t100\t150\tr1\t60\t+\t100\t150\t0,0,0\t1\t50,\t0,\t1\t10\t50\t0.2000000\nchr1\t200\t270\tr2\t30\t-\t200\t270\t0,0,0\t2\t20,20,\t0,50,\t1\t20\t70\t0.2857143\nchr1\t100\t150\tr1\t60\t+\t100\t150\t0,0,0\t1\t50,\t0,\t1\nchr1\t200\t270\tr2\t30\t-\t200\t270\t0,0,0\t2\t20,20,\t0,50,\t0\n' > exp
$BT coverage -a a.tmp.bam -b b.tmp.bed > obs
$BT coverage -a a.tmp.bam -b b.tmp.bed -counts -split >> obs
check obs exp
rm obs exp a.tmp.sam a.tmp.bam b.tmp.bed

###########################################################
#  coverage.new.t020: GFF B is converted from 1-based
#  GFF 101-200 is 0-based 100-200: A 90-110 has 10 of 20
#  bases covered; A 200-210 only touches it.
###########################################################
echo -e "    coverage.new.t020...\c"
printf 'chr1\t90\t110\nchr1\t200\t210\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n' > b.tmp.gff
printf 'chr1\t90\t110\t1\t10\t20\t0.5000000\nchr1\t200\t210\t0\t0\t10\t0.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.gff > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  coverage.new.t021: VCF B spans POS-1 to POS-1+len(REF)
#  On A 120-160, rs1 covers 1 base and rs2 (REF ACGT) 4;
#  rs3 is outside: 2 hits, 5 of 40 bases, 0.1250000.
###########################################################
echo -e "    coverage.new.t021...\c"
printf 'chr1\t120\t160\n' > a.tmp.bed
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t125\trs1\tA\tG\t50\tPASS\tDP=10\nchr1\t150\trs2\tACGT\tA\t50\tPASS\tDP=5\nchr1\t300\trs3\tC\tT\t50\tPASS\tDP=5\n' > b.tmp.vcf
printf 'chr1\t120\t160\t2\t5\t40\t0.1250000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.vcf > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  coverage.new.t022: gzipped A and B read like plain text
###########################################################
echo -e "    coverage.new.t022...\c"
printf 'chr1\t100\t200\tn1\t9\t+\textra1\textra2\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
gzip -c a.tmp.bed > a.tmp.bed.gz
gzip -c b.tmp.bed > b.tmp.bed.gz
printf 'chr1\t100\t200\tn1\t9\t+\textra1\textra2\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed.gz -b b.tmp.bed.gz > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed a.tmp.bed.gz b.tmp.bed.gz

###########################################################
#  coverage.new.t023: -counts with touching, duplicate and
#  zero-length records
#  c1 100-200: B 200-250 only touches; two copies of 150-160
#  overlap -> 2. c2 is zero-length at 300, inside B 280-320
#  -> 1. c3 400-500 holds the zero-length B at 450 -> 1.
#  c4 is zero-length at 600, the point between bases 599 and
#  600; B 598-599 covers only base 598, which is not next to
#  that point -> 0.
###########################################################
echo -e "    coverage.new.t023...\c"
printf 'chr1\t100\t200\tc1\nchr1\t300\t300\tc2\nchr1\t400\t500\tc3\nchr1\t600\t600\tc4\n' > a.tmp.bed
printf 'chr1\t200\t250\nchr1\t150\t160\nchr1\t150\t160\nchr1\t280\t320\nchr1\t450\t450\nchr1\t598\t599\n' > b.tmp.bed
printf 'chr1\t100\t200\tc1\t2\nchr1\t300\t300\tc2\t1\nchr1\t400\t500\tc3\t1\nchr1\t600\t600\tc4\t0\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -counts > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t024: -counts with -f
#  B 100-150 (0.5 of A) and 0-300 (1.0) pass -f 0.5; B
#  150-160 (0.1) does not: count 2.
###########################################################
echo -e "    coverage.new.t024...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t150\nchr1\t150\t160\nchr1\t0\t300\n' > b.tmp.bed
printf 'chr1\t100\t200\t2\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -counts -f 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t025: -d with overlapping B
#  A 100-106; B 101-103, 102-105, 104-110. Depth at
#  100..105 is 0, 1, 2, 1, 2, 1, reported at positions 1-6.
###########################################################
echo -e "    coverage.new.t025...\c"
printf 'chr1\t100\t106\n' > a.tmp.bed
printf 'chr1\t101\t103\nchr1\t102\t105\nchr1\t104\t110\n' > b.tmp.bed
printf 'chr1\t100\t106\t1\t0\nchr1\t100\t106\t2\t1\nchr1\t100\t106\t3\t2\nchr1\t100\t106\t4\t1\nchr1\t100\t106\t5\t2\nchr1\t100\t106\t6\t1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t026: -d on two A records, one uncovered
#  a (chr1:0-3) has base 1 covered; b (chr2:0-2) has none.
###########################################################
echo -e "    coverage.new.t026...\c"
printf 'chr1\t0\t3\ta\nchr2\t0\t2\tb\n' > a.tmp.bed
printf 'chr1\t1\t2\n' > b.tmp.bed
printf 'chr1\t0\t3\ta\t1\t0\nchr1\t0\t3\ta\t2\t1\nchr1\t0\t3\ta\t3\t0\nchr2\t0\t2\tb\t1\t0\nchr2\t0\t2\tb\t2\t0\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t027: -hist per record and the all lines
#  h1 (chr1:0-4) under B 0-2 and 1-3 has depths 1,2,1,0;
#  h2 (chr2:0-2) under B 0-2 has depths 1,1. Per-record lines
#  list each depth in ascending order with its share of A.
#  all: 6 bases, depth 0 x1, 1 x4, 2 x1 -> 0.1666667,
#  0.6666667, 0.1666667.
###########################################################
echo -e "    coverage.new.t027...\c"
printf 'chr1\t0\t4\th1\nchr2\t0\t2\th2\n' > a.tmp.bed
printf 'chr1\t0\t2\nchr1\t1\t3\nchr2\t0\t2\n' > b.tmp.bed
printf 'chr1\t0\t4\th1\t0\t1\t4\t0.2500000\nchr1\t0\t4\th1\t1\t2\t4\t0.5000000\nchr1\t0\t4\th1\t2\t1\t4\t0.2500000\nchr2\t0\t2\th2\t1\t2\t2\t1.0000000\nall\t0\t1\t6\t0.1666667\nall\t1\t4\t6\t0.6666667\nall\t2\t1\t6\t0.1666667\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -hist > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t028: -hist on an A with no coverage
#  All 5 bases at depth 0, in the record line and all line.
###########################################################
echo -e "    coverage.new.t028...\c"
printf 'chr1\t0\t5\n' > a.tmp.bed
printf 'chr2\t0\t5\n' > b.tmp.bed
printf 'chr1\t0\t5\t0\t5\t5\t1.0000000\nall\t0\t5\t5\t1.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -hist > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t029: -hist counts overlapping A separately
#  o1 0-10 and o2 5-15 each have 5 bases under B 5-10. The
#  all lines total 20 bases (10 at depth 0, 10 at depth 1).
###########################################################
echo -e "    coverage.new.t029...\c"
printf 'chr1\t0\t10\to1\nchr1\t5\t15\to2\n' > a.tmp.bed
printf 'chr1\t5\t10\n' > b.tmp.bed
printf 'chr1\t0\t10\to1\t0\t5\t10\t0.5000000\nchr1\t0\t10\to1\t1\t5\t10\t0.5000000\nchr1\t5\t15\to2\t0\t5\t10\t0.5000000\nchr1\t5\t15\to2\t1\t5\t10\t0.5000000\nall\t0\t10\t20\t0.5000000\nall\t1\t10\t20\t0.5000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -hist > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t030: -mean for no, partial and double depth
#  m1 (10 bp) has no B: 0.0000000. m2 (3 bp) has one base at
#  depth 1: 1/3 = 0.3333333. m3 (2 bp) has two copies of B
#  over both bases: 2.0000000.
###########################################################
echo -e "    coverage.new.t030...\c"
printf 'chr1\t0\t10\tm1\nchr1\t20\t23\tm2\nchr1\t30\t32\tm3\n' > a.tmp.bed
printf 'chr1\t22\t23\nchr1\t30\t32\nchr1\t30\t32\n' > b.tmp.bed
printf 'chr1\t0\t10\tm1\t0.0000000\nchr1\t20\t23\tm2\t0.3333333\nchr1\t30\t32\tm3\t2.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -mean > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t031: -mean is rounded as a float
#  B 0-7 and 0-3 on A 0-7 give depth sum 10: 10/7 is
#  1.42857142... in double, but the float nearest it is
#  1.4285714626..., printed as 1.4285715. Pins the current
#  output: float, not double.
###########################################################
echo -e "    coverage.new.t031...\c"
printf 'chr1\t0\t7\n' > a.tmp.bed
printf 'chr1\t0\t7\nchr1\t0\t3\n' > b.tmp.bed
printf 'chr1\t0\t7\t1.4285715\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -mean > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t032: -s counts only same-strand B
#  A on +, - and unknown (.). B: bp 110-120 (+), bm 130-150
#  (-), bd 160-165 (.). p gets bp (10 bases), m gets bm (20),
#  d gets nothing: an unknown strand never matches, as stated
#  in Record::sameChromIntersects.
###########################################################
echo -e "    coverage.new.t032...\c"
printf 'chr1\t100\t200\tp\t0\t+\nchr1\t100\t200\tm\t0\t-\nchr1\t100\t200\td\t0\t.\n' > a.tmp.bed
printf 'chr1\t110\t120\tbp\t0\t+\nchr1\t130\t150\tbm\t0\t-\nchr1\t160\t165\tbd\t0\t.\n' > b.tmp.bed
printf 'chr1\t100\t200\tp\t0\t+\t1\t10\t100\t0.1000000\nchr1\t100\t200\tm\t0\t-\t1\t20\t100\t0.2000000\nchr1\t100\t200\td\t0\t.\t0\t0\t100\t0.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -s > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t033: -S counts only opposite-strand B
#  Same input as t032: p gets bm (20), m gets bp (10), d none.
###########################################################
echo -e "    coverage.new.t033...\c"
printf 'chr1\t100\t200\tp\t0\t+\nchr1\t100\t200\tm\t0\t-\nchr1\t100\t200\td\t0\t.\n' > a.tmp.bed
printf 'chr1\t110\t120\tbp\t0\t+\nchr1\t130\t150\tbm\t0\t-\nchr1\t160\t165\tbd\t0\t.\n' > b.tmp.bed
printf 'chr1\t100\t200\tp\t0\t+\t1\t20\t100\t0.2000000\nchr1\t100\t200\tm\t0\t-\t1\t10\t100\t0.1000000\nchr1\t100\t200\td\t0\t.\t0\t0\t100\t0.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -S > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t034: -s with GFF strands
#  B (0-based 120-130) is on -: the + record gets nothing,
#  the - record 10 of 100 bases.
###########################################################
echo -e "    coverage.new.t034...\c"
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=p\nchr1\tsrc\tgene\t101\t200\t.\t-\t.\tID=m\n' > a.tmp.gff
printf 'chr1\t120\t130\tb\t0\t-\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=p\t0\t0\t100\t0.0000000\nchr1\tsrc\tgene\t101\t200\t.\t-\t.\tID=m\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.gff -b b.tmp.bed -s > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  coverage.new.t035: -s with a strandless (BED5) B is an error
###########################################################
echo -e "    coverage.new.t035...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf 'chr1\t150\t160\tb\t5\n' > b.tmp.bed
printf '***** ERROR: stranded coverage requested, but input file b.tmp.bed does not have strands. *****\nexit 1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -s > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t036: -S with VCF input is an error
###########################################################
echo -e "    coverage.new.t036...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t125\trs1\tA\tG\t50\tPASS\tDP=10\n' > b.tmp.vcf
printf '***** ERROR: stranded coverage not supported for VCF file b.tmp.vcf. *****\nexit 1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.vcf -S > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.vcf

###########################################################
#  coverage.new.t037: -s and -S together are an error
###########################################################
echo -e "    coverage.new.t037...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf 'chr1\t150\t160\tb\t0\t+\n' > b.tmp.bed
printf '***** ERROR: request -s for sameStrand, or -S for diffStrand, not both. *****\nexit 1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -s -S > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t038: -f includes a B at exactly the fraction
#  B 100-150 covers 0.5 of A; B 180-190 covers 0.1.
#  -f 0.5 keeps the first (>=): 1 50 100 0.5000000.
#  -f 0.51 keeps neither: 0 0 100 0.0000000.
###########################################################
echo -e "    coverage.new.t038...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t150\nchr1\t180\t190\n' > b.tmp.bed
printf 'chr1\t100\t200\t1\t50\t100\t0.5000000\nchr1\t100\t200\t0\t0\t100\t0.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -f 0.5 > obs
$BT coverage -a a.tmp.bed -b b.tmp.bed -f 0.51 >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t039: -F is a fraction of B
#  B1 120-130 lies wholly in A; B2 150-400 has 0.2 of itself
#  in A. -F 0.5 keeps only B1: 1 hit, 10 bases.
###########################################################
echo -e "    coverage.new.t039...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t120\t130\nchr1\t150\t400\n' > b.tmp.bed
printf 'chr1\t100\t200\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -F 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t040: -f with and without -r
#  B1 150-400: 0.5 of A, 0.2 of B. B2 90-110: 0.1 of A.
#  B3 175-215: 0.25 of A, 0.625 of B.
#  -f 0.25 -r keeps B3 only: 1 hit, 25 bases.
#  -f 0.25 keeps B1 and B3: 2 hits, union 150-200 = 50 bases.
###########################################################
echo -e "    coverage.new.t040...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t400\nchr1\t90\t110\nchr1\t175\t215\n' > b.tmp.bed
printf 'chr1\t100\t200\t1\t25\t100\t0.2500000\nchr1\t100\t200\t2\t50\t100\t0.5000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -f 0.25 -r > obs
$BT coverage -a a.tmp.bed -b b.tmp.bed -f 0.25 >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t041: -f with -F needs both; -e needs either
#  B1 150-400: 0.5 of A, 0.2 of B1. B2 90-110: 0.1 of A,
#  0.5 of B2.
#  -f 0.5 -F 0.5: neither -> 0 0 100.
#  -f 0.5 -F 0.5 -e: both -> 2 hits, 50 + 10 = 60 bases.
###########################################################
echo -e "    coverage.new.t041...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t400\nchr1\t90\t110\n' > b.tmp.bed
printf 'chr1\t100\t200\t0\t0\t100\t0.0000000\nchr1\t100\t200\t2\t60\t100\t0.6000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.5 > obs
$BT coverage -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.5 -e >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t042: -f and -F outside (0, 1] are errors
#  Checked for -f 0 and -F 1.5.
###########################################################
echo -e "    coverage.new.t042...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: -f must be in the range (0.0, 1.0]. *****\nexit 1\n***** ERROR: -F must be in the range (0.0, 1.0]. *****\nexit 1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -f 0 > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
$BT coverage -a a.tmp.bed -b b.tmp.bed -F 1.5 >> out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t043: -e without -f and -r with -F are errors
###########################################################
echo -e "    coverage.new.t043...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: -e must be used with -f. *****\nexit 1\n***** ERROR: -r must be used solely with -f. *****\nexit 1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -F 0.5 -e > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
$BT coverage -a a.tmp.bed -b b.tmp.bed -F 0.5 -r >> out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t044: BED3 A against a BED12 B
#  B spans 0-300 with blocks 0-110 and 290-300. Without -split
#  all of A 100-200 is covered; with -split only 100-110.
###########################################################
echo -e "    coverage.new.t044...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t0\t300\tbx\t0\t+\t0\t300\t0\t2\t110,10,\t0,290,\n' > b.tmp.bed
printf 'chr1\t100\t200\t1\t100\t100\t1.0000000\nchr1\t100\t200\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
$BT coverage -a a.tmp.bed -b b.tmp.bed -split >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t045: -split with -d and -mean on a BED12 B
#  B covers 0-6 with blocks 0-2 and 4-6: positions 3 and 4 of
#  A 0-6 have depth 0; the mean is 4/6 = 0.6666667.
###########################################################
echo -e "    coverage.new.t045...\c"
printf 'chr1\t0\t6\n' > a.tmp.bed
printf 'chr1\t0\t6\tb\t0\t+\t0\t6\t0\t2\t2,2,\t0,4,\n' > b.tmp.bed
printf 'chr1\t0\t6\t1\t1\nchr1\t0\t6\t2\t1\nchr1\t0\t6\t3\t0\nchr1\t0\t6\t4\t0\nchr1\t0\t6\t5\t1\nchr1\t0\t6\t6\t1\nchr1\t0\t6\t0.6666667\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -split -d > obs
$BT coverage -a a.tmp.bed -b b.tmp.bed -split -mean >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t046: -sorted on several chromosomes
#  chr1 and chr3 have no B; chr4 is only in B.
###########################################################
echo -e "    coverage.new.t046...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\nchr3\t100\t200\n' > a.tmp.bed
printf 'chr2\t150\t160\nchr4\t1\t5\n' > b.tmp.bed
printf 'chr1\t100\t200\t0\t0\t100\t0.0000000\nchr2\t100\t200\t1\t10\t100\t0.1000000\nchr3\t100\t200\t0\t0\t100\t0.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -sorted > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t047: -sorted -g with chr2 ordered before chr1
###########################################################
echo -e "    coverage.new.t047...\c"
printf 'chr2\t100\t200\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr2\t150\t160\nchr1\t120\t140\n' > b.tmp.bed
printf 'chr2\t1000\nchr1\t1000\n' > g.tmp.txt
printf 'chr2\t100\t200\t1\t10\t100\t0.1000000\nchr1\t100\t200\t1\t20\t100\t0.2000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -sorted -g g.tmp.txt > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed g.tmp.txt

###########################################################
#  coverage.new.t048: -sorted with an unsorted A is an error
###########################################################
echo -e "    coverage.new.t048...\c"
printf 'chr1\t300\t400\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'Error: Sorted input specified, but the file a.tmp.bed has the following out of order record\nchr1\t100\t200\nexit 1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -sorted > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t049: -sorted -hist with an A chromosome that
#  B lacks
#  chr1 (4 bp) has no B; chr2 (2 bp) has 1 base covered.
#  all: 6 bases, 5 at depth 0 (0.8333333), 1 at depth 1.
###########################################################
echo -e "    coverage.new.t049...\c"
printf 'chr1\t0\t4\nchr2\t0\t2\n' > a.tmp.bed
printf 'chr2\t0\t1\n' > b.tmp.bed
printf 'chr1\t0\t4\t0\t4\t4\t1.0000000\nchr2\t0\t2\t0\t1\t2\t0.5000000\nchr2\t0\t2\t1\t1\t2\t0.5000000\nall\t0\t5\t6\t0.8333333\nall\t1\t1\t6\t0.1666667\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -sorted -hist > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t050: three -b files, default output
#  b1 110-120, 115-125; b2 150-160, chr2 190-250; b3 118-152.
#  chr1: 4 hits, union 110-160 = 50 bases. chr2: 1 hit, 10.
###########################################################
echo -e "    coverage.new.t050...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t110\t120\nchr1\t115\t125\n' > b1.tmp.bed
printf 'chr1\t150\t160\nchr2\t190\t250\n' > b2.tmp.bed
printf 'chr1\t118\t152\n' > b3.tmp.bed
printf 'chr1\t100\t200\t4\t50\t100\t0.5000000\nchr2\t100\t200\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed b3.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed b3.tmp.bed

###########################################################
#  coverage.new.t051: three -b files with -hist
#  Same input as t050. chr1 depths: 110-114 1, 115-117 2,
#  118-119 3, 120-124 2, 125-149 1, 150-151 2, 152-159 1,
#  else 0 -> depth 0: 50, 1: 38, 2: 10, 3: 2. chr2: 190-199
#  at depth 1. all over 200 bases: 140, 48, 10, 2.
###########################################################
echo -e "    coverage.new.t051...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t110\t120\nchr1\t115\t125\n' > b1.tmp.bed
printf 'chr1\t150\t160\nchr2\t190\t250\n' > b2.tmp.bed
printf 'chr1\t118\t152\n' > b3.tmp.bed
printf 'chr1\t100\t200\t0\t50\t100\t0.5000000\nchr1\t100\t200\t1\t38\t100\t0.3800000\nchr1\t100\t200\t2\t10\t100\t0.1000000\nchr1\t100\t200\t3\t2\t100\t0.0200000\nchr2\t100\t200\t0\t90\t100\t0.9000000\nchr2\t100\t200\t1\t10\t100\t0.1000000\nall\t0\t140\t200\t0.7000000\nall\t1\t48\t200\t0.2400000\nall\t2\t10\t200\t0.0500000\nall\t3\t2\t200\t0.0100000\n' > exp
$BT coverage -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed b3.tmp.bed -hist > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed b3.tmp.bed

###########################################################
#  coverage.new.t052: repeated -b flags act like one list
#  With b1 and b2 of t050: chr1 has 3 hits over 110-125 and
#  150-160 (25 bases); chr2 has 1 hit, 10 bases.
###########################################################
echo -e "    coverage.new.t052...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t110\t120\nchr1\t115\t125\n' > b1.tmp.bed
printf 'chr1\t150\t160\nchr2\t190\t250\n' > b2.tmp.bed
printf 'chr1\t100\t200\t3\t25\t100\t0.2500000\nchr2\t100\t200\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b1.tmp.bed -b b2.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  coverage.new.t053: -d sums depth over several -b files
#  b1 100-103 and b2 102-105 on A 100-105: depths 1,1,2,1,1.
###########################################################
echo -e "    coverage.new.t053...\c"
printf 'chr1\t100\t105\n' > a.tmp.bed
printf 'chr1\t100\t103\n' > b1.tmp.bed
printf 'chr1\t102\t105\n' > b2.tmp.bed
printf 'chr1\t100\t105\t1\t1\nchr1\t100\t105\t2\t1\nchr1\t100\t105\t3\t2\nchr1\t100\t105\t4\t1\nchr1\t100\t105\t5\t1\n' > exp
$BT coverage -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  coverage.new.t054: -header with the default output
#  Without -header the '#' and 'track' lines are dropped; with
#  it they come first, once. B's header is never printed.
###########################################################
echo -e "    coverage.new.t054...\c"
printf '#comment line\ntrack name=foo\nchr1\t100\t200\tr1\nchr2\t10\t20\tr2\n' > a.tmp.bed
printf '#b header\nchr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\tr1\t1\t10\t100\t0.1000000\nchr2\t10\t20\tr2\t0\t0\t10\t0.0000000\n#comment line\ntrack name=foo\nchr1\t100\t200\tr1\t1\t10\t100\t0.1000000\nchr2\t10\t20\tr2\t0\t0\t10\t0.0000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs
$BT coverage -a a.tmp.bed -b b.tmp.bed -header >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t055: -header with -hist
#  The header precedes the record lines and the all lines
#  (110 bases: 100 at depth 0 = 0.9090909, 10 at depth 1 =
#  0.0909091).
###########################################################
echo -e "    coverage.new.t055...\c"
printf '#comment line\ntrack name=foo\nchr1\t100\t200\tr1\nchr2\t10\t20\tr2\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '#comment line\ntrack name=foo\nchr1\t100\t200\tr1\t0\t90\t100\t0.9000000\nchr1\t100\t200\tr1\t1\t10\t100\t0.1000000\nchr2\t10\t20\tr2\t0\t10\t10\t1.0000000\nall\t0\t100\t110\t0.9090909\nall\t1\t10\t110\t0.0909091\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -header -hist > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t056: -header with an A that has only a header
#  The header is printed and the run exits 0.
###########################################################
echo -e "    coverage.new.t056...\c"
printf '#only a header\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '#only a header\nexit 0\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -header > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t057: -counts/-mean and -hist/-mean are
#  mutually exclusive (message and exit status)
###########################################################
echo -e "    coverage.new.t057...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: -counts, -d, -mean, and -hist are all mutually exclusive options. *****\nexit 1\n***** ERROR: -counts, -d, -mean, and -hist are all mutually exclusive options. *****\nexit 1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -counts -mean > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
$BT coverage -a a.tmp.bed -b b.tmp.bed -hist -mean >> out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t058: -a without -b is an error
###########################################################
echo -e "    coverage.new.t058...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf '***** ERROR: query and database files not specified. *****\nexit 1\n' > exp
$BT coverage -a a.tmp.bed > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed

###########################################################
#  coverage.new.t059: a missing input file is an error
#  Checked for a missing A and a missing B.
###########################################################
echo -e "    coverage.new.t059...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'Error: Unable to open file missing.tmp.bed. Exiting.\nexit 1\nError: Unable to open file missing.tmp.bed. Exiting.\nexit 1\n' > exp
$BT coverage -a missing.tmp.bed -b a.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
$BT coverage -a a.tmp.bed -b missing.tmp.bed >> obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  coverage.new.t060: an unknown option is an error
###########################################################
echo -e "    coverage.new.t060...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: Unrecognized parameter: -bogus *****\nexit 1\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -bogus > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t061: chr1 vs 1 naming warns on stderr
#  Nothing overlaps (different names). Pins the current
#  output: B is loaded before A is read, so the warning names
#  A's record, and it is printed when first tripped and again
#  at exit.
###########################################################
echo -e "    coverage.new.t061...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf '1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\t0\t0\t100\t0.0000000\n--- stderr\n***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t100\t200\n\n***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t100\t200\n\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs 2> err.tmp
echo "--- stderr" >> obs
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t062: -nonamecheck silences that warning
###########################################################
echo -e "    coverage.new.t062...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf '1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\t0\t0\t100\t0.0000000\n--- stderr\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -nonamecheck > obs 2> err.tmp
echo "--- stderr" >> obs
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t063: empty A and empty B
#  Empty A prints nothing, even with -hist (no all lines);
#  empty B gives zeros. All runs exit 0.
###########################################################
echo -e "    coverage.new.t063...\c"
printf 'chr1\t100\t200\nchr2\t10\t20\n' > a.tmp.bed
printf '' > e.tmp.bed
printf 'exit 0\nexit 0\nchr1\t100\t200\t0\t0\t100\t0.0000000\nchr2\t10\t20\t0\t0\t10\t0.0000000\nexit 0\n' > exp
$BT coverage -a e.tmp.bed -b a.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
$BT coverage -a e.tmp.bed -b a.tmp.bed -hist >> obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
$BT coverage -a a.tmp.bed -b e.tmp.bed >> obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed e.tmp.bed

###########################################################
#  coverage.new.t064: a VCF structural variant as A, -counts
#  sv2 <DUP> at POS 201 with END=300 covers far more than its
#  1 bp REF. B 202-205 and 250-260 lie inside it under any
#  reading of its first and last base; B 150-190 and 320-330
#  lie well outside: count 2. Its length is not asserted (the
#  SVLEN/END span disagreement is in bugs.md).
###########################################################
echo -e "    coverage.new.t064...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t201\tsv2\tN\t<DUP>\t50\tPASS\tSVTYPE=DUP;END=300\n' > a.tmp.vcf
printf 'chr1\t150\t190\nchr1\t202\t205\nchr1\t250\t260\nchr1\t320\t330\n' > b.tmp.bed
printf 'chr1\t201\tsv2\tN\t<DUP>\t50\tPASS\tSVTYPE=DUP;END=300\t2\n' > exp
$BT coverage -a a.tmp.vcf -b b.tmp.bed -counts > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  coverage.new.t065: -header with BAM A warns and is ignored
#  r1 (100-150) has 10 bases under B 120-130. The warning
#  goes to stderr and no header is printed.
###########################################################
echo -e "    coverage.new.t065...\c"
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*\n' > a.tmp.sam
$htsutil samtobam a.tmp.sam a.tmp.bam
printf 'chr1\t120\t130\n' > b.tmp.bed
printf 'chr1\t100\t150\tr1\t60\t+\t100\t150\t0,0,0\t1\t50,\t0,\t1\t10\t50\t0.2000000\n--- stderr\n\n*****\n*****WARNING: -header option is not valid for BAM input.\n*****\n' > exp
$BT coverage -a a.tmp.bam -b b.tmp.bed -header > obs 2> err.tmp
echo "--- stderr" >> obs
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.sam a.tmp.bam b.tmp.bed

###########################################################
#  coverage.new.t066: -counts accepts BAM A
#  intersect -c refuses BAM A without -bed; coverage always
#  writes text, so -counts works: r1 has 1 overlapping B.
###########################################################
echo -e "    coverage.new.t066...\c"
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*\n' > a.tmp.sam
$htsutil samtobam a.tmp.sam a.tmp.bam
printf 'chr1\t120\t130\n' > b.tmp.bed
printf 'chr1\t100\t150\tr1\t60\t+\t100\t150\t0,0,0\t1\t50,\t0,\t1\n' > exp
$BT coverage -a a.tmp.bam -b b.tmp.bed -counts > obs
check obs exp
rm obs exp a.tmp.sam a.tmp.bam b.tmp.bed

###########################################################
#  coverage.new.t067: -header keeps browser, track and # lines
###########################################################
echo -e "    coverage.new.t067...\c"
printf 'browser position chr1:1-1000\ntrack name=x\n#c\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'browser position chr1:1-1000\ntrack name=x\n#c\nchr1\t100\t200\t1\t10\t100\t0.1000000\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  coverage.new.t068: chr1 vs chr01 naming warns on stderr
#  Nothing overlaps (different names). Pins the current
#  output: B is loaded first, so A's record is named, and the
#  warning is printed when tripped and again at exit.
###########################################################
echo -e "    coverage.new.t068...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr01\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\t0\t0\t100\t0.0000000\n--- stderr\n***** WARNING: File a.tmp.bed has a record where naming convention (leading zero) is inconsistent with other files:\nchr1\t100\t200\n\n***** WARNING: File a.tmp.bed has a record where naming convention (leading zero) is inconsistent with other files:\nchr1\t100\t200\n\n' > exp
$BT coverage -a a.tmp.bed -b b.tmp.bed > obs 2> err.tmp
echo "--- stderr" >> obs
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

[[ $FAILURES -eq 0 ]] || exit 1;
