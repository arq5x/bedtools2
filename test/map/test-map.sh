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
#  Test defaults
############################################################
echo "    map.t01...\c"
echo \
"chr1	0	100	30
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	6
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test sum
############################################################
echo "    map.t02...\c"
echo \
"chr1	0	100	30
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	6
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed -o sum > obs
check obs exp
rm obs exp


###########################################################
#  Test count
############################################################
echo "    map.t03...\c"
echo \
"chr1	0	100	3
chr1	100	200	1
chr2	0	100	0
chr2	100	200	0
chr3	0	100	3
chr3	100	200	1" > exp
$BT map -a ivls.bed -b values.bed -o count > obs
check obs exp
rm obs exp


###########################################################
#  Test mean
############################################################
echo "    map.t04...\c"
echo \
"chr1	0	100	10
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	2
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed -o mean > obs
check obs exp
rm obs exp

###########################################################
#  Test max
############################################################
echo "    map.t05...\c"
echo \
"chr1	0	100	15
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	3
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed -o max > obs
check obs exp
rm obs exp

###########################################################
#  Test min
############################################################
echo "    map.t06...\c"
echo \
"chr1	0	100	5
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	1
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed -o min > obs
check obs exp
rm obs exp

###########################################################
#  Test mode
############################################################
echo "    map.t07...\c"
echo \
"chr1	0	100	5
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	1
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values2.bed -o mode > obs
check obs exp
rm obs exp

###########################################################
#  Test anti-mode
############################################################
echo "    map.t08...\c"
echo \
"chr1	0	100	10
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	1
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values2.bed -o antimode > obs
check obs exp
rm obs exp

