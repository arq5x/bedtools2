set -e;
BT=${BT-../../bin/bedtools}

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
#  Test a basic self intersection
###########################################################
echo -e "    jaccard.t01...\c"
echo \
"intersection	union	jaccard	n_intersections
110	110	1	2" > exp
$BT jaccard -a a.bed -b a.bed > obs
check obs exp
rm obs exp


echo -e "    jaccard.t02...\c"
echo \
"intersection	union	jaccard	n_intersections
10	140	0.0714286	1" > exp
$BT jaccard -a a.bed -b b.bed > obs
check obs exp
rm obs exp

echo -e "    jaccard.t03...\c"
echo \
"intersection	union	jaccard	n_intersections
10	200	0.05	1" > exp
$BT jaccard -a a.bed -b c.bed > obs
check obs exp
rm obs exp

# TEST #4 IS DEPRECATED

###########################################################
#  Test stdin
###########################################################
echo -e "    jaccard.t05...\c"
echo \
"intersection	union	jaccard	n_intersections
10	140	0.0714286	1" > exp
cat a.bed | $BT jaccard -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test symmetry
###########################################################
echo -e "    jaccard.t06...\c"
$BT jaccard -a a.bed -b b.bed > obs1
$BT jaccard -a b.bed -b a.bed > obs2
check obs1 obs2
rm obs1 obs2

###########################################################
#  Test partially matching blocks without -split option.
###########################################################
echo -e "    jaccard.t07...\c"
echo \
"intersection	union	jaccard	n_intersections
10	50	0.2	1" > exp
$BT jaccard -a three_blocks_match.bed -b e.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test partially matching blocks with -split option.
###########################################################
echo -e "    jaccard.t08...\c"
echo \
"intersection	union	jaccard	n_intersections
5	35	0.142857	1" > exp
$BT jaccard -a three_blocks_match.bed -b e.bed -split > obs
check obs exp
rm obs exp

###########################################################
#  Test jaccard of Bam with Bam
###########################################################
echo -e "    jaccard.t09...\c"
echo \
"intersection	union	jaccard	n_intersections
10	150	0.0666667	1" > exp
$BT jaccard -a a.bam -b three_blocks_match.bam -bed > obs
check exp obs
rm exp obs

###########################################################
#  Test jaccard with mixed strand files
###########################################################
echo -e "    jaccard.t10...\c"
echo \
"intersection	union	jaccard	n_intersections
145	180	0.805556	2" >exp
$BT jaccard -a aMixedStrands.bed -b bMixedStrands.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test jaccard with mixed strand files, -s option
#  (match strand, either forward or reverse)
###########################################################
echo -e "    jaccard.t11...\c"
echo \
"intersection	union	jaccard	n_intersections
120	290	0.413793	4" >exp
$BT jaccard -a aMixedStrands.bed -b bMixedStrands.bed -s > obs
check obs exp
rm obs exp

###########################################################
#  Test jaccard with mixed strand files, -S + option
#  (match strand, forward only)
###########################################################
echo -e "    jaccard.t12...\c"
echo \
"intersection	union	jaccard	n_intersections
40	135	0.296296	2" >exp
$BT jaccard -a aMixedStrands.bed -b bMixedStrands.bed -S + > obs
check obs exp
rm obs exp

###########################################################
#  Test jaccard with mixed strand files, -S - option
#  (match strand, reverse only)
###########################################################
echo -e "    jaccard.t13...\c"
echo \
"intersection	union	jaccard	n_intersections
80	155	0.516129	2" > exp
$BT jaccard -a aMixedStrands.bed -b bMixedStrands.bed -S - > obs
check obs exp
rm obs exp

echo -e "    jaccard.t14...\c"
echo  "intersection	union	jaccard	n_intersections
1	3	0.333333	1" > exp

$BT jaccard -b a645.bed -a b645.bed > obs
check obs exp
rm obs 

echo -e "    jaccard.t15...\c"
$BT jaccard -a a645.bed -b b645.bed > obs
check obs exp
rm obs  exp

echo -e "    jaccard.t16...\c"
echo -e "intersection	union	jaccard	n_intersections
247800000	615800000	0.402403	4" > exp
$BT jaccard -a long.bed -b short.bed > obs
check obs exp
rm obs  exp

[[ $FAILURES -eq 0 ]] || exit 1;
