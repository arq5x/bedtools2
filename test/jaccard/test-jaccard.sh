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
#  Test a basic self intersection
###########################################################
echo "    jaccard.t01...\c"
echo \
"intersection	union	jaccard	n_intersections
110	110	1	2" > exp
$BT jaccard -a a.bed -b a.bed > obs
check obs exp
rm obs exp


echo "    jaccard.t02...\c"
echo \
"intersection	union	jaccard	n_intersections
10	140	0.0714286	1" > exp
$BT jaccard -a a.bed -b b.bed > obs
check obs exp
rm obs exp

echo "    jaccard.t03...\c"
echo \
"intersection	union	jaccard	n_intersections
10	200	0.05	1" > exp
$BT jaccard -a a.bed -b c.bed > obs
check obs exp
rm obs exp


echo "    jaccard.t04...\c"
echo \
"intersection	union	jaccard	n_intersections
0	210	0	0" > exp
$BT jaccard -a a.bed -b d.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test stdin
###########################################################
echo "    jaccard.t05...\c"
echo \
"intersection	union	jaccard	n_intersections
10	140	0.0714286	1" > exp
cat a.bed | $BT jaccard -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test symmetry
###########################################################
echo "    jaccard.t06...\c"
$BT jaccard -a a.bed -b b.bed > obs1
$BT jaccard -a b.bed -b a.bed > obs2
check obs1 obs2

