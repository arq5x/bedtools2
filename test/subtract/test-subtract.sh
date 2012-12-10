BT=${BT-../../bin/bedtools}

check()
{
	if diff $1 $2; then
    	echo ok
		return 1
	else
    	echo fail
		return 0
	fi
}

# cat a.bed
# chr1	10	20	a1	1	+
# chr1	50	70	a2	2	-
# 
# cat b.bed
# chr1	18	25	b1	1	-
# chr1	80	90	b2	2	+


###########################################################
# test baseline subtraction
###########################################################
echo "    subtract.t1...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo "    subtract.t2...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.5 > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo "    subtract.t3...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.1 > obs
check obs exp
rm obs exp


###########################################################
# test -s subtraction
###########################################################
echo "    subtract.t4...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -s > obs
check obs exp
rm obs exp


###########################################################
# test -S subtraction
###########################################################
echo "    subtract.t5...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -S > obs
check obs exp
rm obs exp


###########################################################
# test -A subtraction
###########################################################
echo "    subtract.t6...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A > obs
check obs exp
rm obs exp


###########################################################
# test -A with -f subtraction
###########################################################
echo "    subtract.t7...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.5 > obs
check obs exp
rm obs exp

###########################################################
# test -A with -f subtraction
###########################################################
echo "    subtract.t8...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.1 > obs
check obs exp
rm obs exp


###########################################################
# test -N with -f subtraction
###########################################################
echo "    subtract.t9...\c"
echo \
"chr1	0	10" > exp
$BT subtract -a c.bed -b d.bed -N -f 0.4 > obs
check obs exp
rm obs exp

###########################################################
# test -N with -f subtraction
###########################################################
echo "    subtract.t10...\c"
touch exp
$BT subtract -a c.bed -b d.bed -N -f 0.39 > obs
check obs exp
rm obs exp