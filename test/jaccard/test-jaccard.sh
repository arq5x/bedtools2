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
#  Test a basic self intersection
############################################################
echo "    jaccard.t01...\c"
echo \
"intersection	union	jaccard
110	110	1" > exp
$BT jaccard -a a.bed -b a.bed > obs
check obs exp
rm obs exp


echo "    jaccard.t02...\c"
echo \
"intersection	union	jaccard
10	140	0.0714286" > exp
$BT jaccard -a a.bed -b b.bed > obs
check obs exp
rm obs exp

echo "    jaccard.t03...\c"
echo \
"intersection	union	jaccard
10	200	0.05" > exp
$BT jaccard -a a.bed -b c.bed > obs
check obs exp
rm obs exp


echo "    jaccard.t04...\c"
echo \
"intersection	union	jaccard
0	210	0" > exp
$BT jaccard -a a.bed -b b.bed > obs
check obs exp
rm obs exp