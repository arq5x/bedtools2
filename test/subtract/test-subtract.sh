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
echo -e "    subtract.t1...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo -e "    subtract.t2...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.5 > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo -e "    subtract.t3...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.1 > obs
check obs exp
rm obs exp


###########################################################
# test -s subtraction
###########################################################
echo -e "    subtract.t4...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -s > obs
check obs exp
rm obs exp


###########################################################
# test -S subtraction
###########################################################
echo -e "    subtract.t5...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -S > obs
check obs exp
rm obs exp


###########################################################
# test -A subtraction
###########################################################
echo -e "    subtract.t6...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A > obs
check obs exp
rm obs exp


###########################################################
# test -A with -f subtraction
###########################################################
echo -e "    subtract.t7...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.5 > obs
check obs exp
rm obs exp

###########################################################
# test -A with -f subtraction
###########################################################
echo -e "    subtract.t8...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.1 > obs
check obs exp
rm obs exp


###########################################################
# test -N with -f subtraction
###########################################################
echo -e "    subtract.t9...\c"
echo \
"chr1	0	10" > exp
$BT subtract -a c.bed -b d.bed -N -f 0.4 > obs
check obs exp
rm obs exp

###########################################################
# test -N with -f subtraction
###########################################################
echo -e "    subtract.t10...\c"
touch exp
$BT subtract -a c.bed -b d.bed -N -f 0.39 > obs
check obs exp
rm obs exp


###########################################################
# 
# REPEAT TESTS 1-8 WITH -sorted 
#
###########################################################


###########################################################
# test baseline subtraction
###########################################################
echo -e "    subtract.t11...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo -e "    subtract.t12...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.5  -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo -e "    subtract.t13...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.1  -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -s subtraction
###########################################################
echo -e "    subtract.t14...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -s -sorted  > obs
check obs exp
rm obs exp


###########################################################
# test -S subtraction
###########################################################
echo -e "    subtract.t15...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -S -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -A subtraction
###########################################################
echo -e "    subtract.t16...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A  -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -A with -f subtraction
###########################################################
echo -e "    subtract.t17...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.5  -sorted > obs
check obs exp
rm obs exp

###########################################################
# test -A with -f subtraction
###########################################################
echo -e "    subtract.t18...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.1  -sorted > obs
check obs exp
rm obs exp


###########################################################
#
# TEST WITH MULTIPLE DATABASES
# 
###########################################################


###########################################################
# test with 2 DBs
###########################################################
echo -e "    subtract.t19...\c"
echo \
"chr1	15	18	a1	1	+
chr1	50	55	a2	2	-
chr1	65	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed b2.bed > obs
check obs exp
rm obs exp

###########################################################
# test with 2 DBs, -f option
###########################################################
echo -e "    subtract.t20...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed b2.bed -f 0.8 > obs
check obs exp
rm obs exp

###########################################################
# test with 2 DBs, -f option and -N options
###########################################################
echo -e "    subtract.t21...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed b2.bed -f 0.6 -N > obs
check obs exp
rm obs exp

[[ $FAILURES -eq 0 ]] || exit 1;
