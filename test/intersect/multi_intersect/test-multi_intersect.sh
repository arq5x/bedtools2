set -e;
echo -e \
"\n###########################################################
#  
#  MULTIPLE DATABASE INTERSECTION
#
###########################################################\n"


BT=${BT-../../../bin/bedtools}

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
#  Test the intersection of a query with 3 dbs against 
#  the merged and sorted result of 3 
#  seperate intersections of the query with each db.
#
#  First make sure the 3 seperate intersections are correct.
#
############################################################
echo -e "    intersect.t01...\c"
echo \
"chr1	5	20
chr1	70	75
chr2	5	20
chr2	70	75
chr3	5	20
chr3	70	75" > exp1 
$BT intersect -a query.bed -b d1.bed > obs
check obs exp1
rm obs

echo -e "    intersect.t02...\c"
echo \
"chr1	40	45
chr1	110	120
chr2	40	45
chr2	110	120
chr3	40	45
chr3	110	120" > exp2
$BT intersect -a query.bed -b d2.bed > obs
check obs exp2
rm obs

echo -e "    intersect.t03...\c"
echo \
"chr1	85	90
chr1	105	115
chr2	85	90
chr2	105	115
chr3	85	90
chr3	105	115" > exp3
$BT intersect -a query.bed -b d3.bed > obs
check obs exp3
rm obs

###########################################################
#  Now test that the R-tree (unsorted input) with sorted
#  output will give the correct result
###########################################################
echo -e "    intersect.t04...\c"
cat exp1 exp2 exp3 | sort -k1,1 -k2,2n > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -sortout > obs
check exp obs
rm obs exp1 exp2 exp3


###########################################################
#  And then test that sorted input with sorted
#  output will give the correct result
###########################################################
echo -e "    intersect.t05...\c"
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -sorted -sortout > obs
check exp obs
rm obs exp


###########################################################
# 
#  Now repeat the above with chroms in a non-standard order
#
###########################################################

echo -e "    intersect.t06...\c"
echo \
"chr1	5	20
chr1	70	75
chr2	5	20
chr2	70	75
chr10	5	20
chr10	70	75" > exp1 
$BT intersect -a query2.bed -b d4.bed > obs
check obs exp1
rm obs

echo -e "    intersect.t07...\c"
echo \
"chr1	40	45
chr1	110	120
chr2	40	45
chr2	110	120
chr10	40	45
chr10	110	120" > exp2
$BT intersect -a query2.bed -b d5.bed > obs
check obs exp2
rm obs

echo -e "    intersect.t08...\c"
echo \
"chr1	85	90
chr1	105	115
chr2	85	90
chr2	105	115
chr10	85	90
chr10	105	115" > exp3
$BT intersect -a query2.bed -b d6.bed > obs
check obs exp3
rm obs

###########################################################
#  Now test that the R-tree (unsorted input) with sorted
#  output will give the correct result
###########################################################
echo -e "    intersect.t09...\c"
cat exp1 exp2 exp3 | sort -k1.4,1n -k2,2n > exp
$BT intersect -a query2.bed -b d4.bed d5.bed d6.bed -g g.bed -sortout > obs
check exp obs
rm obs exp1 exp2 exp3


###########################################################
#  And then test that sorted input with sorted
#  output will give the correct result
###########################################################
echo -e "    intersect.t10...\c"
$BT intersect -a query2.bed -b d4.bed d5.bed d6.bed -sorted -g g.bed -sortout > obs
check exp obs
rm obs exp



###########################################################
# 
# TEST VARIOUS CMD LINE OUTPUT OPTIONS
#
###########################################################


###########################################################
#  Test -c option
###########################################################
echo -e "    intersect.t11...\c"
echo \
"chr1	1	20	1
chr1	40	45	1
chr1	70	90	2
chr1	105	120	2
chr2	1	20	1
chr2	40	45	1
chr2	70	90	2
chr2	105	120	2
chr3	1	20	1
chr3	40	45	1
chr3	70	90	2
chr3	105	120	2" > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -c -sorted -sortout > obs
check exp obs
rm exp obs

###########################################################
#  Test -C option
###########################################################
echo -e "    intersect.t12...\c"
echo \
"chr1	1	20	1	1
chr1	1	20	2	0
chr1	1	20	3	0
chr1	40	45	1	0
chr1	40	45	2	1
chr1	40	45	3	0
chr1	70	90	1	1
chr1	70	90	2	0
chr1	70	90	3	1
chr1	105	120	1	0
chr1	105	120	2	1
chr1	105	120	3	1
chr2	1	20	1	1
chr2	1	20	2	0
chr2	1	20	3	0
chr2	40	45	1	0
chr2	40	45	2	1
chr2	40	45	3	0
chr2	70	90	1	1
chr2	70	90	2	0
chr2	70	90	3	1
chr2	105	120	1	0
chr2	105	120	2	1
chr2	105	120	3	1
chr3	1	20	1	1
chr3	1	20	2	0
chr3	1	20	3	0
chr3	40	45	1	0
chr3	40	45	2	1
chr3	40	45	3	0
chr3	70	90	1	1
chr3	70	90	2	0
chr3	70	90	3	1
chr3	105	120	1	0
chr3	105	120	2	1
chr3	105	120	3	1" > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -C > obs
check exp obs
rm exp obs

###########################################################
#  Test -C with -filenames option
###########################################################
echo -e "    intersect.t13...\c"
echo \
"chr1	1	20	d1.bed	1
chr1	1	20	d2.bed	0
chr1	1	20	d3.bed	0
chr1	40	45	d1.bed	0
chr1	40	45	d2.bed	1
chr1	40	45	d3.bed	0
chr1	70	90	d1.bed	1
chr1	70	90	d2.bed	0
chr1	70	90	d3.bed	1
chr1	105	120	d1.bed	0
chr1	105	120	d2.bed	1
chr1	105	120	d3.bed	1
chr2	1	20	d1.bed	1
chr2	1	20	d2.bed	0
chr2	1	20	d3.bed	0
chr2	40	45	d1.bed	0
chr2	40	45	d2.bed	1
chr2	40	45	d3.bed	0
chr2	70	90	d1.bed	1
chr2	70	90	d2.bed	0
chr2	70	90	d3.bed	1
chr2	105	120	d1.bed	0
chr2	105	120	d2.bed	1
chr2	105	120	d3.bed	1
chr3	1	20	d1.bed	1
chr3	1	20	d2.bed	0
chr3	1	20	d3.bed	0
chr3	40	45	d1.bed	0
chr3	40	45	d2.bed	1
chr3	40	45	d3.bed	0
chr3	70	90	d1.bed	1
chr3	70	90	d2.bed	0
chr3	70	90	d3.bed	1
chr3	105	120	d1.bed	0
chr3	105	120	d2.bed	1
chr3	105	120	d3.bed	1" > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -C -filenames > obs
check exp obs
rm exp obs

###########################################################
#  Test -C with -names option
###########################################################
echo -e "    intersect.t14...\c"
echo \
"chr1	1	20	d1	1
chr1	1	20	d2	0
chr1	1	20	d3	0
chr1	40	45	d1	0
chr1	40	45	d2	1
chr1	40	45	d3	0
chr1	70	90	d1	1
chr1	70	90	d2	0
chr1	70	90	d3	1
chr1	105	120	d1	0
chr1	105	120	d2	1
chr1	105	120	d3	1
chr2	1	20	d1	1
chr2	1	20	d2	0
chr2	1	20	d3	0
chr2	40	45	d1	0
chr2	40	45	d2	1
chr2	40	45	d3	0
chr2	70	90	d1	1
chr2	70	90	d2	0
chr2	70	90	d3	1
chr2	105	120	d1	0
chr2	105	120	d2	1
chr2	105	120	d3	1
chr3	1	20	d1	1
chr3	1	20	d2	0
chr3	1	20	d3	0
chr3	40	45	d1	0
chr3	40	45	d2	1
chr3	40	45	d3	0
chr3	70	90	d1	1
chr3	70	90	d2	0
chr3	70	90	d3	1
chr3	105	120	d1	0
chr3	105	120	d2	1
chr3	105	120	d3	1" > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -C -names d1 d2 d3 > obs
check exp obs
rm exp obs

###########################################################
#  Test -C with -names option and -f 0.5
###########################################################
echo -e "    intersect.t15...\c"
echo \
"chr1	1	20	d1	1
chr1	1	20	d2	0
chr1	1	20	d3	0
chr1	40	45	d1	0
chr1	40	45	d2	1
chr1	40	45	d3	0
chr1	70	90	d1	0
chr1	70	90	d2	0
chr1	70	90	d3	0
chr1	105	120	d1	0
chr1	105	120	d2	1
chr1	105	120	d3	1
chr2	1	20	d1	1
chr2	1	20	d2	0
chr2	1	20	d3	0
chr2	40	45	d1	0
chr2	40	45	d2	1
chr2	40	45	d3	0
chr2	70	90	d1	0
chr2	70	90	d2	0
chr2	70	90	d3	0
chr2	105	120	d1	0
chr2	105	120	d2	1
chr2	105	120	d3	1
chr3	1	20	d1	1
chr3	1	20	d2	0
chr3	1	20	d3	0
chr3	40	45	d1	0
chr3	40	45	d2	1
chr3	40	45	d3	0
chr3	70	90	d1	0
chr3	70	90	d2	0
chr3	70	90	d3	0
chr3	105	120	d1	0
chr3	105	120	d2	1
chr3	105	120	d3	1" > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -C -f 0.5 -names d1 d2 d3 > obs
check exp obs
rm exp obs

###########################################################
#  Test -C with -names option and -f 1.0
###########################################################
echo -e "    intersect.t16...\c"
echo \
"chr1	1	20	d1	0
chr1	1	20	d2	0
chr1	1	20	d3	0
chr1	40	45	d1	0
chr1	40	45	d2	1
chr1	40	45	d3	0
chr1	70	90	d1	0
chr1	70	90	d2	0
chr1	70	90	d3	0
chr1	105	120	d1	0
chr1	105	120	d2	0
chr1	105	120	d3	0
chr2	1	20	d1	0
chr2	1	20	d2	0
chr2	1	20	d3	0
chr2	40	45	d1	0
chr2	40	45	d2	1
chr2	40	45	d3	0
chr2	70	90	d1	0
chr2	70	90	d2	0
chr2	70	90	d3	0
chr2	105	120	d1	0
chr2	105	120	d2	0
chr2	105	120	d3	0
chr3	1	20	d1	0
chr3	1	20	d2	0
chr3	1	20	d3	0
chr3	40	45	d1	0
chr3	40	45	d2	1
chr3	40	45	d3	0
chr3	70	90	d1	0
chr3	70	90	d2	0
chr3	70	90	d3	0
chr3	105	120	d1	0
chr3	105	120	d2	0
chr3	105	120	d3	0" > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -C -f 1.0 -names d1 d2 d3 > obs
check exp obs
rm exp obs


###########################################################
#  Test -f option
###########################################################
echo -e "    intersect.t17...\c"
echo \
"chr1	5	20
chr1	40	45
chr1	110	120
chr1	105	115
chr2	5	20
chr2	40	45
chr2	110	120
chr2	105	115
chr3	5	20
chr3	40	45
chr3	110	120
chr3	105	115" > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -f .6 > obs
check exp obs


###########################################################
#  Test -s option
###########################################################
echo -e "    intersect.t18...\c"
echo \
"chr1	5	20	q1	100	+
chr1	85	90	q3	100	+
chr1	110	120	q4	100	+
chr1	105	115	q4	100	+
chr3	5	20	q9	100	+
chr3	85	90	q11	100	+
chr3	110	120	q12	100	+
chr3	105	115	q12	100	+" > exp
$BT intersect -a query3.bed -b d7.bed d8.bed d9.bed -s > obs
check exp obs
rm exp obs


###########################################################
#  Test -s option
###########################################################
echo -e "    intersect.t13...\c"
echo \
"chr1	40	45	q2	100	+
chr1	70	75	q3	100	+
chr2	5	20	q5	100	+
chr2	40	45	q6	100	+
chr2	85	90	q7	100	+
chr2	105	115	q8	100	+
chr3	40	45	q10	100	+
chr3	70	75	q11	100	+" > exp
$BT intersect -a query3.bed -b d7.bed d8.bed d9.bed -S > obs
check exp obs
rm exp obs


###########################################################
#  Test -wo option
###########################################################
echo -e "    intersect.t19...\c"
echo \
"chr1	1	20	1	chr1	5	25	15
chr1	40	45	2	chr1	40	50	5
chr1	70	90	1	chr1	65	75	5
chr1	70	90	3	chr1	85	115	5
chr1	105	120	2	chr1	110	125	10
chr1	105	120	3	chr1	85	115	10
chr2	1	20	1	chr2	5	25	15
chr2	40	45	2	chr2	40	50	5
chr2	70	90	1	chr2	65	75	5
chr2	70	90	3	chr2	85	115	5
chr2	105	120	2	chr2	110	125	10
chr2	105	120	3	chr2	85	115	10
chr3	1	20	1	chr3	5	25	15
chr3	40	45	2	chr3	40	50	5
chr3	70	90	1	chr3	65	75	5
chr3	70	90	3	chr3	85	115	5
chr3	105	120	2	chr3	110	125	10
chr3	105	120	3	chr3	85	115	10" > exp
$BT intersect -a query.bed -b d1.bed d2.bed d3.bed -wo > obs
check exp obs
rm exp obs



###########################################################
#  Test -wa -wb -header option
###########################################################
echo -e "    intersect.t20...\c"
echo \
"# Query 3 Header Line.
chr1	1	20	q1	100	+	1	chr1	5	25	d7_1	100	+
chr1	40	45	q2	100	+	2	chr1	40	50	d8_1	100	-
chr1	70	90	q3	100	+	1	chr1	65	75	d7_2	100	-
chr1	70	90	q3	100	+	3	chr1	85	115	d9_1	100	+
chr1	105	120	q4	100	+	2	chr1	110	125	d8_2	100	+
chr1	105	120	q4	100	+	3	chr1	85	115	d9_1	100	+
chr2	1	20	q5	100	+	1	chr2	5	25	d7_4	100	-
chr2	40	45	q6	100	+	2	chr2	40	50	d8_3	100	-
chr2	70	90	q7	100	+	1	chr2	65	75	d7_5	100	.
chr2	70	90	q7	100	+	3	chr2	85	115	d9_1	100	-
chr2	105	120	q8	100	+	2	chr2	110	125	d8_4	100	.
chr2	105	120	q8	100	+	3	chr2	85	115	d9_1	100	-
chr3	1	20	q9	100	+	1	chr3	5	25	d7_7	100	+
chr3	40	45	q10	100	+	2	chr3	40	50	d8_5	100	-
chr3	70	90	q11	100	+	1	chr3	65	75	d7_8	100	-
chr3	70	90	q11	100	+	3	chr3	85	115	d9_1	100	+
chr3	105	120	q12	100	+	2	chr3	110	125	d8_6	100	+
chr3	105	120	q12	100	+	3	chr3	85	115	d9_1	100	+" > exp
$BT intersect -a query3.bed -b d7.bed d8.bed d9.bed -wa -wb -header > obs
check exp obs
rm exp obs



###########################################################
#  Test the -filenames option, before db listing
###########################################################
echo -e "    intersect.t21...\c"
echo \
"chr1	1	20	q1	100	+	d7.bed	chr1	5	25	d7_1	100	+
chr1	40	45	q2	100	+	d8.bed	chr1	40	50	d8_1	100	-
chr1	70	90	q3	100	+	d7.bed	chr1	65	75	d7_2	100	-
chr1	70	90	q3	100	+	d9.bed	chr1	85	115	d9_1	100	+
chr1	105	120	q4	100	+	d8.bed	chr1	110	125	d8_2	100	+
chr1	105	120	q4	100	+	d9.bed	chr1	85	115	d9_1	100	+
chr2	1	20	q5	100	+	d7.bed	chr2	5	25	d7_4	100	-
chr2	40	45	q6	100	+	d8.bed	chr2	40	50	d8_3	100	-
chr2	70	90	q7	100	+	d7.bed	chr2	65	75	d7_5	100	.
chr2	70	90	q7	100	+	d9.bed	chr2	85	115	d9_1	100	-
chr2	105	120	q8	100	+	d8.bed	chr2	110	125	d8_4	100	.
chr2	105	120	q8	100	+	d9.bed	chr2	85	115	d9_1	100	-
chr3	1	20	q9	100	+	d7.bed	chr3	5	25	d7_7	100	+
chr3	40	45	q10	100	+	d8.bed	chr3	40	50	d8_5	100	-
chr3	70	90	q11	100	+	d7.bed	chr3	65	75	d7_8	100	-
chr3	70	90	q11	100	+	d9.bed	chr3	85	115	d9_1	100	+
chr3	105	120	q12	100	+	d8.bed	chr3	110	125	d8_6	100	+
chr3	105	120	q12	100	+	d9.bed	chr3	85	115	d9_1	100	+" > exp
$BT intersect -a query3.bed -filenames -b d7.bed d8.bed d9.bed -wa -wb > obs
check exp obs
rm exp obs

###########################################################
#  Test the -filenames option, after db listing
###########################################################
echo -e "    intersect.t22...\c"
echo \
"chr1	1	20	q1	100	+	d7.bed	chr1	5	25	d7_1	100	+
chr1	40	45	q2	100	+	d8.bed	chr1	40	50	d8_1	100	-
chr1	70	90	q3	100	+	d7.bed	chr1	65	75	d7_2	100	-
chr1	70	90	q3	100	+	d9.bed	chr1	85	115	d9_1	100	+
chr1	105	120	q4	100	+	d8.bed	chr1	110	125	d8_2	100	+
chr1	105	120	q4	100	+	d9.bed	chr1	85	115	d9_1	100	+
chr2	1	20	q5	100	+	d7.bed	chr2	5	25	d7_4	100	-
chr2	40	45	q6	100	+	d8.bed	chr2	40	50	d8_3	100	-
chr2	70	90	q7	100	+	d7.bed	chr2	65	75	d7_5	100	.
chr2	70	90	q7	100	+	d9.bed	chr2	85	115	d9_1	100	-
chr2	105	120	q8	100	+	d8.bed	chr2	110	125	d8_4	100	.
chr2	105	120	q8	100	+	d9.bed	chr2	85	115	d9_1	100	-
chr3	1	20	q9	100	+	d7.bed	chr3	5	25	d7_7	100	+
chr3	40	45	q10	100	+	d8.bed	chr3	40	50	d8_5	100	-
chr3	70	90	q11	100	+	d7.bed	chr3	65	75	d7_8	100	-
chr3	70	90	q11	100	+	d9.bed	chr3	85	115	d9_1	100	+
chr3	105	120	q12	100	+	d8.bed	chr3	110	125	d8_6	100	+
chr3	105	120	q12	100	+	d9.bed	chr3	85	115	d9_1	100	+" > exp
$BT intersect -a query3.bed -b d7.bed d8.bed d9.bed -filenames -wa -wb > obs
check exp obs
rm exp obs


###########################################################
#  Test the -names option, before db listing
###########################################################
echo -e "    intersect.t23...\c"
echo \
"chr1	1	20	q1	100	+	blue	chr1	5	25	d7_1	100	+
chr1	40	45	q2	100	+	red	chr1	40	50	d8_1	100	-
chr1	70	90	q3	100	+	blue	chr1	65	75	d7_2	100	-
chr1	70	90	q3	100	+	green	chr1	85	115	d9_1	100	+
chr1	105	120	q4	100	+	red	chr1	110	125	d8_2	100	+
chr1	105	120	q4	100	+	green	chr1	85	115	d9_1	100	+
chr2	1	20	q5	100	+	blue	chr2	5	25	d7_4	100	-
chr2	40	45	q6	100	+	red	chr2	40	50	d8_3	100	-
chr2	70	90	q7	100	+	blue	chr2	65	75	d7_5	100	.
chr2	70	90	q7	100	+	green	chr2	85	115	d9_1	100	-
chr2	105	120	q8	100	+	red	chr2	110	125	d8_4	100	.
chr2	105	120	q8	100	+	green	chr2	85	115	d9_1	100	-
chr3	1	20	q9	100	+	blue	chr3	5	25	d7_7	100	+
chr3	40	45	q10	100	+	red	chr3	40	50	d8_5	100	-
chr3	70	90	q11	100	+	blue	chr3	65	75	d7_8	100	-
chr3	70	90	q11	100	+	green	chr3	85	115	d9_1	100	+
chr3	105	120	q12	100	+	red	chr3	110	125	d8_6	100	+
chr3	105	120	q12	100	+	green	chr3	85	115	d9_1	100	+" > exp
$BT intersect -a query3.bed -names blue red green -b d7.bed d8.bed d9.bed -wa -wb > obs
check exp obs
rm exp obs

###########################################################
#  Test the -names option, after db listing
###########################################################
echo -e "    intersect.t24...\c"
echo \
"chr1	1	20	q1	100	+	blue	chr1	5	25	d7_1	100	+
chr1	40	45	q2	100	+	red	chr1	40	50	d8_1	100	-
chr1	70	90	q3	100	+	blue	chr1	65	75	d7_2	100	-
chr1	70	90	q3	100	+	green	chr1	85	115	d9_1	100	+
chr1	105	120	q4	100	+	red	chr1	110	125	d8_2	100	+
chr1	105	120	q4	100	+	green	chr1	85	115	d9_1	100	+
chr2	1	20	q5	100	+	blue	chr2	5	25	d7_4	100	-
chr2	40	45	q6	100	+	red	chr2	40	50	d8_3	100	-
chr2	70	90	q7	100	+	blue	chr2	65	75	d7_5	100	.
chr2	70	90	q7	100	+	green	chr2	85	115	d9_1	100	-
chr2	105	120	q8	100	+	red	chr2	110	125	d8_4	100	.
chr2	105	120	q8	100	+	green	chr2	85	115	d9_1	100	-
chr3	1	20	q9	100	+	blue	chr3	5	25	d7_7	100	+
chr3	40	45	q10	100	+	red	chr3	40	50	d8_5	100	-
chr3	70	90	q11	100	+	blue	chr3	65	75	d7_8	100	-
chr3	70	90	q11	100	+	green	chr3	85	115	d9_1	100	+
chr3	105	120	q12	100	+	red	chr3	110	125	d8_6	100	+
chr3	105	120	q12	100	+	green	chr3	85	115	d9_1	100	+" > exp
$BT intersect -a query3.bed -b d7.bed d8.bed d9.bed -names blue red green -wa -wb > obs
check exp obs
rm exp obs



[[ $FAILURES -eq 0 ]] || exit 1;
