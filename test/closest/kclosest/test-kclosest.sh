echo -e \
"\n###########################################################
#  
#  K CLOSEST HITS TESTS
#
###########################################################\n"

BT=${BT-../../../bin/bedtools}

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

###########################################################
# Test that hits are correctly picked for various values
# of k and tie modes.
###########################################################
echo "    kclosest.t1...\c"
echo \
"chr1	100	110	chr1	70	80
chr1	100	110	chr1	50	60
chr1	100	110	chr1	30	40
chr1	100	110	chr1	35	40
chr1	100	110	chr1	38	40" > exp
$BT closest -a q1.bed -b d1.bed -k 3 > obs
check obs exp
rm obs exp


echo "    kclosest.t2...\c"
echo \
"chr1	100	110	chr1	70	80
chr1	100	110	chr1	50	60
chr1	100	110	chr1	30	40
chr1	100	110	chr1	35	40
chr1	100	110	chr1	38	40" > exp
$BT closest -a q1.bed -b d1.bed -k 5 > obs
check obs exp
rm obs exp


echo "    kclosest.t3...\c"
echo \
"chr1	100	110	chr1	70	80
chr1	100	110	chr1	50	60
chr1	100	110	chr1	30	40
chr1	100	110	chr1	35	40
chr1	100	110	chr1	38	40
chr1	100	110	chr1	25	28" > exp
$BT closest -a q1.bed -b d1.bed -k 6 > obs
check obs exp
rm obs exp


echo "    kclosest.t4...\c"
echo \
"chr1	100	110	chr1	70	80
chr1	100	110	chr1	50	60
chr1	100	110	chr1	30	40
chr1	100	110	chr1	35	40
chr1	100	110	chr1	38	40
chr1	100	110	chr1	25	28
chr1	100	110	chr1	10	15" > exp
$BT closest -a q1.bed -b d1.bed -k 7 > obs
check obs exp
rm obs exp


echo "    kclosest.t5...\c"
echo \
"chr1	100	110	chr1	70	80
chr1	100	110	chr1	50	60
chr1	100	110	chr1	30	40
chr1	100	110	chr1	25	28" > exp
$BT closest -a q1.bed -b d1.bed -k 4 -t first > obs
check obs exp
rm obs exp


echo "    kclosest.t6...\c"
echo \
"chr1	100	110	chr1	70	80
chr1	100	110	chr1	50	60
chr1	100	110	chr1	38	40
chr1	100	110	chr1	25	28" > exp
$BT closest -a q1.bed -b d1.bed -k 4 -t last > obs
check obs exp
rm obs exp


echo "    kclosest.t7...\c"
echo \
"chr1	100	110	chr1	95	105	2.1	20	+
chr1	100	110	chr1	98	108	2.2	20	+
chr1	100	110	chr1	105	115	2.3	20	+
chr1	100	110	chr1	120	130	2.4	20	+" > exp
$BT closest -a q1.bed -b d2.bed -k 4 > obs
check obs exp
rm obs exp

echo "    kclosest.t8...\c"
echo \
"chr1	100	110	chr1	95	105	2.1	20	+
chr1	100	110	chr1	120	130	2.4	20	+
chr1	100	110	chr1	140	160	2.5	20	+
chr1	100	110	chr1	170	180	2.8	20	+" > exp
$BT closest -a q1.bed -b d2.bed -k 4 -t first > obs
check obs exp
rm obs exp


echo "    kclosest.t9...\c"
echo \
"chr1	100	110	chr1	120	130	2.4	20	+
chr1	100	110	chr1	140	160	2.5	20	+
chr1	100	110	chr1	170	180	2.8	20	+
chr1	100	110	chr1	190	200	2.9	20	+" > exp
$BT closest -a q1.bed -b d2.bed -k 4 -t first -io > obs
check obs exp
rm obs exp

echo "    kclosest.t10...\c"
echo \
"chr1	100	110	chr1	95	105	2.1	20	+
chr1	100	110	chr1	98	108	2.2	20	+
chr1	100	110	chr1	105	115	2.3	20	+
chr1	100	110	chr1	120	130	2.4	20	+
chr1	100	110	chr1	140	160	2.5	20	+
chr1	100	110	chr1	140	155	2.6	20	+
chr1	100	110	chr1	140	150	2.7	20	+" > exp
$BT closest -a q1.bed -b d2.bed -k 7  > obs
check obs exp
rm obs exp

echo "    kclosest.t11...\c"
echo \
"chr1	100	110	chr1	95	105	2.1	20	+
chr1	100	110	chr1	120	130	2.4	20	+
chr1	100	110	chr1	140	160	2.5	20	+
chr1	100	110	chr1	170	180	2.8	20	+
chr1	100	110	chr1	190	200	2.9	20	+" > exp
$BT closest -a q1.bed -b d2.bed -k 7 -t first > obs
check obs exp
rm obs exp

echo "    kclosest.t12...\c"
echo \
"chr1	100	110	chr1	105	115	2.3	20	+
chr1	100	110	chr1	120	130	2.4	20	+
chr1	100	110	chr1	140	150	2.7	20	+
chr1	100	110	chr1	170	180	2.8	20	+
chr1	100	110	chr1	190	200	2.9	20	+" > exp
$BT closest -a q1.bed -b d2.bed -k 7 -t last > obs
check obs exp
rm obs exp


echo "    kclosest.t13...\c"
echo \
"chr1	100	110	chr1	95	105	2.1	20	+	0
chr1	100	110	chr1	98	108	2.2	20	+	0
chr1	100	110	chr1	105	115	2.3	20	+	0" > exp
$BT closest -a q1.bed -b d2.bed -k 7 -id -D ref > obs
check obs exp
rm obs exp

echo "    kclosest.t14...\c"
echo \
"chr1	100	110	chr1	95	105	2.1	20	+	0" > exp
$BT closest -a q1.bed -b d2.bed -k 7 -id -D ref -t first > obs
check obs exp
rm obs exp


###########################################################
# Check strandedness with -iu, -io, -id, and -D ref,
# a, and b.
###########################################################
echo "    kclosest.t15...\c"
echo \
"chr1	100	110	chr1	90	105	3.7	20	+	0
chr1	100	110	chr1	95	110	3.8	20	-	0
chr1	100	110	chr1	110	115	3.9	20	+	1
chr1	100	110	chr1	120	130	3.10	20	-	11
chr1	100	110	chr1	130	140	3.11	20	+	21
chr1	100	110	chr1	135	140	3.115	20	+	26
chr1	100	110	chr1	150	165	3.12	20	-	41
chr1	100	110	chr1	150	160	3.13	20	-	41
chr1	100	110	chr1	150	155	3.14	20	-	41
chr1	100	110	chr1	170	180	3.15	20	+	61
chr1	100	110	chr1	190	200	3.16	20	-	81" > exp
$BT closest -a q1.bed -b d3.bed -k 15 -D ref -iu > obs
check obs exp
rm obs exp

echo "    kclosest.t16...\c"
echo \
"chr1	100	110	chr1	110	115	3.9	20	+	1
chr1	100	110	chr1	120	130	3.10	20	-	11
chr1	100	110	chr1	70	80	3.6	20	+	-21
chr1	100	110	chr1	75	80	3.65	20	+	-21
chr1	100	110	chr1	130	140	3.11	20	+	21
chr1	100	110	chr1	135	140	3.115	20	+	26
chr1	100	110	chr1	45	60	3.3	20	-	-41
chr1	100	110	chr1	50	60	3.4	20	-	-41
chr1	100	110	chr1	55	60	3.5	20	-	-41
chr1	100	110	chr1	150	165	3.12	20	-	41
chr1	100	110	chr1	150	160	3.13	20	-	41
chr1	100	110	chr1	150	155	3.14	20	-	41
chr1	100	110	chr1	30	40	3.2	20	+	-61
chr1	100	110	chr1	170	180	3.15	20	+	61
chr1	100	110	chr1	10	20	3.1	20	-	-81
chr1	100	110	chr1	190	200	3.16	20	-	81" > exp
$BT closest -a q1.bed -b d3.bed -k 15 -D ref -io > obs
check obs exp
rm obs exp


echo "    kclosest.t17...\c"
echo \
"chr1	100	110	chr1	90	105	3.7	20	+	0
chr1	100	110	chr1	95	110	3.8	20	-	0
chr1	100	110	chr1	70	80	3.6	20	+	-21
chr1	100	110	chr1	75	80	3.65	20	+	-21
chr1	100	110	chr1	45	60	3.3	20	-	-41
chr1	100	110	chr1	50	60	3.4	20	-	-41
chr1	100	110	chr1	55	60	3.5	20	-	-41
chr1	100	110	chr1	30	40	3.2	20	+	-61
chr1	100	110	chr1	10	20	3.1	20	-	-81" > exp
$BT closest -a q1.bed -b d3.bed -k 15 -D ref -id > obs
check obs exp
rm obs exp

echo "    kclosest.t18...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	chr1	110	115	3.9	20	+	1
chr1	100	110	q2.1	20	+	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	chr1	130	140	3.11	20	+	21
chr1	100	110	q2.1	20	+	chr1	135	140	3.115	20	+	26
chr1	100	110	q2.1	20	+	chr1	150	165	3.12	20	-	41
chr1	100	110	q2.1	20	+	chr1	150	160	3.13	20	-	41
chr1	100	110	q2.1	20	+	chr1	150	155	3.14	20	-	41
chr1	100	110	q2.1	20	+	chr1	170	180	3.15	20	+	61
chr1	100	110	q2.1	20	+	chr1	190	200	3.16	20	-	81
chr1	105	110	q2.2	20	-	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	chr1	90	105	3.7	20	+	1
chr1	105	110	q2.2	20	-	chr1	70	80	3.6	20	+	26
chr1	105	110	q2.2	20	-	chr1	75	80	3.65	20	+	26
chr1	105	110	q2.2	20	-	chr1	45	60	3.3	20	-	46
chr1	105	110	q2.2	20	-	chr1	50	60	3.4	20	-	46
chr1	105	110	q2.2	20	-	chr1	55	60	3.5	20	-	46
chr1	105	110	q2.2	20	-	chr1	30	40	3.2	20	+	66
chr1	105	110	q2.2	20	-	chr1	10	20	3.1	20	-	86" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D a -iu > obs
check obs exp
rm obs exp

echo "    kclosest.t19...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	chr1	70	80	3.6	20	+	-21
chr1	100	110	q2.1	20	+	chr1	75	80	3.65	20	+	-21
chr1	100	110	q2.1	20	+	chr1	45	60	3.3	20	-	-41
chr1	100	110	q2.1	20	+	chr1	50	60	3.4	20	-	-41
chr1	100	110	q2.1	20	+	chr1	55	60	3.5	20	-	-41
chr1	100	110	q2.1	20	+	chr1	30	40	3.2	20	+	-61
chr1	100	110	q2.1	20	+	chr1	10	20	3.1	20	-	-81
chr1	105	110	q2.2	20	-	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	chr1	110	115	3.9	20	+	-1
chr1	105	110	q2.2	20	-	chr1	120	130	3.10	20	-	-11
chr1	105	110	q2.2	20	-	chr1	130	140	3.11	20	+	-21
chr1	105	110	q2.2	20	-	chr1	135	140	3.115	20	+	-26
chr1	105	110	q2.2	20	-	chr1	150	165	3.12	20	-	-41
chr1	105	110	q2.2	20	-	chr1	150	160	3.13	20	-	-41
chr1	105	110	q2.2	20	-	chr1	150	155	3.14	20	-	-41
chr1	105	110	q2.2	20	-	chr1	170	180	3.15	20	+	-61
chr1	105	110	q2.2	20	-	chr1	190	200	3.16	20	-	-81" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D a -id > obs
check obs exp
rm obs exp


echo "    kclosest.t20...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	chr1	70	80	3.6	20	+	21
chr1	100	110	q2.1	20	+	chr1	75	80	3.65	20	+	21
chr1	100	110	q2.1	20	+	chr1	150	165	3.12	20	-	41
chr1	100	110	q2.1	20	+	chr1	150	160	3.13	20	-	41
chr1	100	110	q2.1	20	+	chr1	150	155	3.14	20	-	41
chr1	100	110	q2.1	20	+	chr1	30	40	3.2	20	+	61
chr1	100	110	q2.1	20	+	chr1	190	200	3.16	20	-	81
chr1	105	110	q2.2	20	-	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	chr1	90	105	3.7	20	+	1
chr1	105	110	q2.2	20	-	chr1	120	130	3.10	20	-	11
chr1	105	110	q2.2	20	-	chr1	70	80	3.6	20	+	26
chr1	105	110	q2.2	20	-	chr1	75	80	3.65	20	+	26
chr1	105	110	q2.2	20	-	chr1	150	165	3.12	20	-	41
chr1	105	110	q2.2	20	-	chr1	150	160	3.13	20	-	41
chr1	105	110	q2.2	20	-	chr1	150	155	3.14	20	-	41
chr1	105	110	q2.2	20	-	chr1	30	40	3.2	20	+	66
chr1	105	110	q2.2	20	-	chr1	190	200	3.16	20	-	81" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D b -iu > obs
check obs exp
rm obs exp

echo "    kclosest.t21...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	chr1	110	115	3.9	20	+	-1
chr1	100	110	q2.1	20	+	chr1	130	140	3.11	20	+	-21
chr1	100	110	q2.1	20	+	chr1	135	140	3.115	20	+	-26
chr1	100	110	q2.1	20	+	chr1	45	60	3.3	20	-	-41
chr1	100	110	q2.1	20	+	chr1	50	60	3.4	20	-	-41
chr1	100	110	q2.1	20	+	chr1	55	60	3.5	20	-	-41
chr1	100	110	q2.1	20	+	chr1	170	180	3.15	20	+	-61
chr1	100	110	q2.1	20	+	chr1	10	20	3.1	20	-	-81
chr1	105	110	q2.2	20	-	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	chr1	110	115	3.9	20	+	-1
chr1	105	110	q2.2	20	-	chr1	130	140	3.11	20	+	-21
chr1	105	110	q2.2	20	-	chr1	135	140	3.115	20	+	-26
chr1	105	110	q2.2	20	-	chr1	45	60	3.3	20	-	-46
chr1	105	110	q2.2	20	-	chr1	50	60	3.4	20	-	-46
chr1	105	110	q2.2	20	-	chr1	55	60	3.5	20	-	-46
chr1	105	110	q2.2	20	-	chr1	170	180	3.15	20	+	-61
chr1	105	110	q2.2	20	-	chr1	10	20	3.1	20	-	-86" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D b -id > obs
check obs exp
rm obs exp

echo "    kclosest.t22...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	chr1	110	115	3.9	20	+	1
chr1	100	110	q2.1	20	+	chr1	130	140	3.11	20	+	21
chr1	100	110	q2.1	20	+	chr1	135	140	3.115	20	+	26
chr1	100	110	q2.1	20	+	chr1	170	180	3.15	20	+	61
chr1	105	110	q2.2	20	-	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	chr1	45	60	3.3	20	-	46
chr1	105	110	q2.2	20	-	chr1	50	60	3.4	20	-	46
chr1	105	110	q2.2	20	-	chr1	55	60	3.5	20	-	46
chr1	105	110	q2.2	20	-	chr1	10	20	3.1	20	-	86" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D a -s -iu > obs
check obs exp
rm obs exp

echo "    kclosest.t23...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	chr1	70	80	3.6	20	+	-21
chr1	100	110	q2.1	20	+	chr1	75	80	3.65	20	+	-21
chr1	100	110	q2.1	20	+	chr1	30	40	3.2	20	+	-61
chr1	105	110	q2.2	20	-	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	chr1	120	130	3.10	20	-	-11
chr1	105	110	q2.2	20	-	chr1	150	165	3.12	20	-	-41
chr1	105	110	q2.2	20	-	chr1	150	160	3.13	20	-	-41
chr1	105	110	q2.2	20	-	chr1	150	155	3.14	20	-	-41
chr1	105	110	q2.2	20	-	chr1	190	200	3.16	20	-	-81" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D a -s -id > obs
check obs exp
rm obs exp

echo "    kclosest.t24...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	chr1	70	80	3.6	20	+	21
chr1	100	110	q2.1	20	+	chr1	75	80	3.65	20	+	21
chr1	100	110	q2.1	20	+	chr1	30	40	3.2	20	+	61
chr1	105	110	q2.2	20	-	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	chr1	120	130	3.10	20	-	11
chr1	105	110	q2.2	20	-	chr1	150	165	3.12	20	-	41
chr1	105	110	q2.2	20	-	chr1	150	160	3.13	20	-	41
chr1	105	110	q2.2	20	-	chr1	150	155	3.14	20	-	41
chr1	105	110	q2.2	20	-	chr1	190	200	3.16	20	-	81" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D b -s -iu > obs
check obs exp
rm obs exp

echo "    kclosest.t25...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	chr1	110	115	3.9	20	+	-1
chr1	100	110	q2.1	20	+	chr1	130	140	3.11	20	+	-21
chr1	100	110	q2.1	20	+	chr1	135	140	3.115	20	+	-26
chr1	100	110	q2.1	20	+	chr1	170	180	3.15	20	+	-61
chr1	105	110	q2.2	20	-	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	chr1	45	60	3.3	20	-	-46
chr1	105	110	q2.2	20	-	chr1	50	60	3.4	20	-	-46
chr1	105	110	q2.2	20	-	chr1	55	60	3.5	20	-	-46
chr1	105	110	q2.2	20	-	chr1	10	20	3.1	20	-	-86" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D b -s -id > obs
check obs exp
rm obs exp

echo "    kclosest.t26...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	chr1	150	165	3.12	20	-	41
chr1	100	110	q2.1	20	+	chr1	150	160	3.13	20	-	41
chr1	100	110	q2.1	20	+	chr1	150	155	3.14	20	-	41
chr1	100	110	q2.1	20	+	chr1	190	200	3.16	20	-	81
chr1	105	110	q2.2	20	-	chr1	90	105	3.7	20	+	1
chr1	105	110	q2.2	20	-	chr1	70	80	3.6	20	+	26
chr1	105	110	q2.2	20	-	chr1	75	80	3.65	20	+	26
chr1	105	110	q2.2	20	-	chr1	30	40	3.2	20	+	66" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D a -S -iu > obs
check obs exp
rm obs exp

echo "    kclosest.t27...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	chr1	45	60	3.3	20	-	-41
chr1	100	110	q2.1	20	+	chr1	50	60	3.4	20	-	-41
chr1	100	110	q2.1	20	+	chr1	55	60	3.5	20	-	-41
chr1	100	110	q2.1	20	+	chr1	10	20	3.1	20	-	-81
chr1	105	110	q2.2	20	-	chr1	110	115	3.9	20	+	-1
chr1	105	110	q2.2	20	-	chr1	130	140	3.11	20	+	-21
chr1	105	110	q2.2	20	-	chr1	135	140	3.115	20	+	-26
chr1	105	110	q2.2	20	-	chr1	170	180	3.15	20	+	-61" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D a -S -id > obs
check obs exp
rm obs exp

echo "    kclosest.t28...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	chr1	150	165	3.12	20	-	41
chr1	100	110	q2.1	20	+	chr1	150	160	3.13	20	-	41
chr1	100	110	q2.1	20	+	chr1	150	155	3.14	20	-	41
chr1	100	110	q2.1	20	+	chr1	190	200	3.16	20	-	81
chr1	105	110	q2.2	20	-	chr1	90	105	3.7	20	+	1
chr1	105	110	q2.2	20	-	chr1	70	80	3.6	20	+	26
chr1	105	110	q2.2	20	-	chr1	75	80	3.65	20	+	26
chr1	105	110	q2.2	20	-	chr1	30	40	3.2	20	+	66" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D b -S -iu > obs
check obs exp
rm obs exp

echo "    kclosest.t29...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	chr1	45	60	3.3	20	-	-41
chr1	100	110	q2.1	20	+	chr1	50	60	3.4	20	-	-41
chr1	100	110	q2.1	20	+	chr1	55	60	3.5	20	-	-41
chr1	100	110	q2.1	20	+	chr1	10	20	3.1	20	-	-81
chr1	105	110	q2.2	20	-	chr1	110	115	3.9	20	+	-1
chr1	105	110	q2.2	20	-	chr1	130	140	3.11	20	+	-21
chr1	105	110	q2.2	20	-	chr1	135	140	3.115	20	+	-26
chr1	105	110	q2.2	20	-	chr1	170	180	3.15	20	+	-61" > exp
$BT closest -a q2.bed -b d3.bed -k 15 -D b -S -id > obs
check obs exp
rm obs exp

###########################################################
# Test multiple databases and mdb = all mode
###########################################################


echo "    kclosest.t30...\c"
echo \
"chr1	100	110	q2.1	20	+	2	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	1	chr1	95	105	2.1	20	+	0
chr1	100	110	q2.1	20	+	2	chr1	95	110	3.8	20	-	0
chr1	100	110	q2.1	20	+	1	chr1	98	108	2.2	20	+	0
chr1	100	110	q2.1	20	+	1	chr1	105	115	2.3	20	+	0
chr1	100	110	q2.1	20	+	2	chr1	110	115	3.9	20	+	1
chr1	100	110	q2.1	20	+	1	chr1	120	130	2.4	20	+	11
chr1	100	110	q2.1	20	+	2	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	2	chr1	70	80	3.6	20	+	-21
chr1	100	110	q2.1	20	+	2	chr1	75	80	3.65	20	+	-21
chr1	100	110	q2.1	20	+	2	chr1	130	140	3.11	20	+	21
chr1	105	110	q2.2	20	-	2	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	1	chr1	98	108	2.2	20	+	0
chr1	105	110	q2.2	20	-	1	chr1	105	115	2.3	20	+	0
chr1	105	110	q2.2	20	-	2	chr1	90	105	3.7	20	+	-1
chr1	105	110	q2.2	20	-	1	chr1	95	105	2.1	20	+	-1
chr1	105	110	q2.2	20	-	2	chr1	110	115	3.9	20	+	1
chr1	105	110	q2.2	20	-	1	chr1	120	130	2.4	20	+	11
chr1	105	110	q2.2	20	-	2	chr1	120	130	3.10	20	-	11
chr1	105	110	q2.2	20	-	2	chr1	130	140	3.11	20	+	21
chr1	105	110	q2.2	20	-	2	chr1	70	80	3.6	20	+	-26
chr1	105	110	q2.2	20	-	2	chr1	75	80	3.65	20	+	-26
chr1	105	110	q2.2	20	-	2	chr1	135	140	3.115	20	+	26" > exp
$BT closest -a q2.bed -b d2.bed d3.bed -mdb all -k 10  -D ref -t all > obs
check obs exp
rm obs exp

echo "    kclosest.t31...\c"
echo \
"chr1	100	110	q2.1	20	+	2	chr1	90	105	3.7	20	+	0
chr1	100	110	q2.1	20	+	2	chr1	110	115	3.9	20	+	1
chr1	100	110	q2.1	20	+	1	chr1	120	130	2.4	20	+	11
chr1	100	110	q2.1	20	+	2	chr1	70	80	3.6	20	+	-21
chr1	100	110	q2.1	20	+	2	chr1	135	140	3.115	20	+	26
chr1	100	110	q2.1	20	+	1	chr1	140	160	2.5	20	+	31
chr1	100	110	q2.1	20	+	2	chr1	45	60	3.3	20	-	-41
chr1	100	110	q2.1	20	+	2	chr1	30	40	3.2	20	+	-61
chr1	100	110	q2.1	20	+	2	chr1	10	20	3.1	20	-	-81
chr1	105	110	q2.2	20	-	2	chr1	95	110	3.8	20	-	0
chr1	105	110	q2.2	20	-	2	chr1	90	105	3.7	20	+	-1
chr1	105	110	q2.2	20	-	1	chr1	120	130	2.4	20	+	11
chr1	105	110	q2.2	20	-	2	chr1	130	140	3.11	20	+	21
chr1	105	110	q2.2	20	-	2	chr1	70	80	3.6	20	+	-26
chr1	105	110	q2.2	20	-	1	chr1	140	160	2.5	20	+	31
chr1	105	110	q2.2	20	-	2	chr1	150	165	3.12	20	-	41
chr1	105	110	q2.2	20	-	2	chr1	45	60	3.3	20	-	-46
chr1	105	110	q2.2	20	-	1	chr1	170	180	2.8	20	+	61
chr1	105	110	q2.2	20	-	2	chr1	30	40	3.2	20	+	-66" > exp
$BT closest -a q2.bed -b d2.bed d3.bed -mdb all -k 10  -D ref -t first > obs
check obs exp
rm obs exp


echo "    kclosest.t32...\c"
echo \
"chr1	100	110	q2.1	20	+	1	chr1	105	115	2.3	20	+	0
chr1	100	110	q2.1	20	+	2	chr1	110	115	3.9	20	+	1
chr1	100	110	q2.1	20	+	2	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	2	chr1	130	140	3.11	20	+	21
chr1	100	110	q2.1	20	+	2	chr1	135	140	3.115	20	+	26
chr1	100	110	q2.1	20	+	1	chr1	140	150	2.7	20	+	31
chr1	100	110	q2.1	20	+	2	chr1	150	155	3.14	20	-	41
chr1	100	110	q2.1	20	+	2	chr1	170	180	3.15	20	+	61
chr1	100	110	q2.1	20	+	2	chr1	190	200	3.16	20	-	81
chr1	105	110	q2.2	20	-	1	chr1	105	115	2.3	20	+	0
chr1	105	110	q2.2	20	-	2	chr1	110	115	3.9	20	+	1
chr1	105	110	q2.2	20	-	2	chr1	120	130	3.10	20	-	11
chr1	105	110	q2.2	20	-	2	chr1	130	140	3.11	20	+	21
chr1	105	110	q2.2	20	-	2	chr1	135	140	3.115	20	+	26
chr1	105	110	q2.2	20	-	1	chr1	140	150	2.7	20	+	31
chr1	105	110	q2.2	20	-	2	chr1	150	155	3.14	20	-	41
chr1	105	110	q2.2	20	-	2	chr1	55	60	3.5	20	-	-46
chr1	105	110	q2.2	20	-	2	chr1	170	180	3.15	20	+	61
chr1	105	110	q2.2	20	-	2	chr1	30	40	3.2	20	+	-66" > exp
$BT closest -a q2.bed -b d2.bed d3.bed -mdb all -k 10  -D ref -t last > obs
check obs exp
rm obs exp

###########################################################
# Test new -fu and -fd features
###########################################################

echo "    kclosest.t33...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	70	80	3.6	20	+	-21
chr1	100	110	q2.1	20	+	chr1	45	60	3.3	20	-	-41
chr1	100	110	q2.1	20	+	chr1	30	40	3.2	20	+	-61
chr1	100	110	q2.1	20	+	chr1	10	20	3.1	20	-	-81
chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	105	110	q2.2	20	-	chr1	90	105	3.7	20	+	-1
chr1	105	110	q2.2	20	-	chr1	70	80	3.6	20	+	-26
chr1	105	110	q2.2	20	-	chr1	45	60	3.3	20	-	-46
chr1	105	110	q2.2	20	-	chr1	30	40	3.2	20	+	-66
chr1	105	110	q2.2	20	-	chr1	10	20	3.1	20	-	-86" > exp
$BT closest -a q2.bed -b d3.bed -k 5  -D ref -t first -fu > obs
check obs exp
rm obs exp

echo "    kclosest.t34...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	110	115	3.9	20	+	1
chr1	100	110	q2.1	20	+	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	chr1	130	140	3.11	20	+	21
chr1	100	110	q2.1	20	+	chr1	135	140	3.115	20	+	26
chr1	100	110	q2.1	20	+	chr1	150	165	3.12	20	-	41
chr1	105	110	q2.2	20	-	chr1	110	115	3.9	20	+	1
chr1	105	110	q2.2	20	-	chr1	120	130	3.10	20	-	11
chr1	105	110	q2.2	20	-	chr1	130	140	3.11	20	+	21
chr1	105	110	q2.2	20	-	chr1	135	140	3.115	20	+	26
chr1	105	110	q2.2	20	-	chr1	150	165	3.12	20	-	41" > exp
$BT closest -a q2.bed -b d3.bed -k 5  -D ref -t first -fd > obs
check obs exp
rm obs exp

echo "    kclosest.t35...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	70	80	3.6	20	+	-21
chr1	100	110	q2.1	20	+	chr1	45	60	3.3	20	-	-41
chr1	100	110	q2.1	20	+	chr1	30	40	3.2	20	+	-61
chr1	100	110	q2.1	20	+	chr1	10	20	3.1	20	-	-81
chr1	100	110	q2.1	20	+	chr1	90	105	3.7	20	+	0
chr1	105	110	q2.2	20	-	chr1	110	115	3.9	20	+	-1
chr1	105	110	q2.2	20	-	chr1	120	130	3.10	20	-	-11
chr1	105	110	q2.2	20	-	chr1	130	140	3.11	20	+	-21
chr1	105	110	q2.2	20	-	chr1	135	140	3.115	20	+	-26
chr1	105	110	q2.2	20	-	chr1	150	165	3.12	20	-	-41" > exp
$BT closest -a q2.bed -b d3.bed -k 5  -D a -t first -fu > obs
check obs exp
rm obs exp

echo "    kclosest.t36...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	110	115	3.9	20	+	1
chr1	100	110	q2.1	20	+	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	chr1	130	140	3.11	20	+	21
chr1	100	110	q2.1	20	+	chr1	135	140	3.115	20	+	26
chr1	100	110	q2.1	20	+	chr1	150	165	3.12	20	-	41
chr1	105	110	q2.2	20	-	chr1	90	105	3.7	20	+	1
chr1	105	110	q2.2	20	-	chr1	70	80	3.6	20	+	26
chr1	105	110	q2.2	20	-	chr1	45	60	3.3	20	-	46
chr1	105	110	q2.2	20	-	chr1	30	40	3.2	20	+	66
chr1	105	110	q2.2	20	-	chr1	10	20	3.1	20	-	86" > exp
$BT closest -a q2.bed -b d3.bed -k 5  -D a -t first -fd > obs
check obs exp
rm obs exp

echo "    kclosest.t37...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	110	115	3.9	20	+	-1
chr1	100	110	q2.1	20	+	chr1	130	140	3.11	20	+	-21
chr1	100	110	q2.1	20	+	chr1	135	140	3.115	20	+	-26
chr1	100	110	q2.1	20	+	chr1	45	60	3.3	20	-	-41
chr1	100	110	q2.1	20	+	chr1	170	180	3.15	20	+	-61
chr1	105	110	q2.2	20	-	chr1	110	115	3.9	20	+	-1
chr1	105	110	q2.2	20	-	chr1	130	140	3.11	20	+	-21
chr1	105	110	q2.2	20	-	chr1	135	140	3.115	20	+	-26
chr1	105	110	q2.2	20	-	chr1	45	60	3.3	20	-	-46
chr1	105	110	q2.2	20	-	chr1	170	180	3.15	20	+	-61" > exp
$BT closest -a q2.bed -b d3.bed -k 5  -D b -t first -fu > obs
check obs exp
rm obs exp

echo "    kclosest.t38...\c"
echo \
"chr1	100	110	q2.1	20	+	chr1	120	130	3.10	20	-	11
chr1	100	110	q2.1	20	+	chr1	70	80	3.6	20	+	21
chr1	100	110	q2.1	20	+	chr1	150	165	3.12	20	-	41
chr1	100	110	q2.1	20	+	chr1	30	40	3.2	20	+	61
chr1	100	110	q2.1	20	+	chr1	190	200	3.16	20	-	81
chr1	105	110	q2.2	20	-	chr1	90	105	3.7	20	+	1
chr1	105	110	q2.2	20	-	chr1	120	130	3.10	20	-	11
chr1	105	110	q2.2	20	-	chr1	70	80	3.6	20	+	26
chr1	105	110	q2.2	20	-	chr1	150	165	3.12	20	-	41
chr1	105	110	q2.2	20	-	chr1	30	40	3.2	20	+	66" > exp
$BT closest -a q2.bed -b d3.bed -k 5  -D b -t first -fd > obs
check obs exp
rm obs exp



