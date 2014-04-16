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
# chr1	10	20
# chr1	30	40
# chr1	40	50
# chr1	45	100

###########################################################
#  Test a basic merge; one interval should be un-merged, 
#  the other two should be merged.
###########################################################
echo "    merge.t1...\c"
echo \
"chr1	10	20
chr1	30	100" > exp
$BT merge -i a.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test that -n option is shown as deperecated
###########################################################
echo "    merge.t2...\c"
echo "***** ERROR: -n option is deprecated. Please see the documentation for the -c and -o column operation options. *****" > exp
$BT merge -i a.bed -n 2>&1 > /dev/null | head -2 | tail -1 > obs
check obs exp
rm obs exp


###########################################################
#  Test the counting of merged intervals. (-n)
###########################################################
echo "    merge.t3...\c"
echo \
"chr1	10	20	1
chr1	30	100	3" > exp
$BT merge -i a.bed -c 1 -o count > obs
check obs exp
rm obs exp


###########################################################
#  Test that -nms option is deprecated
###########################################################
echo "    merge.t4...\c"
echo "***** ERROR: -nms option is deprecated. Please see the documentation for the -c and -o column operation options. *****" > exp
$BT merge -i a.bed -nms 2>&1 > /dev/null | head -2 | tail -1 > obs
check obs exp
rm obs exp

###########################################################
#  Test the listing of names from merged intervals.
###########################################################
echo "    merge.t5...\c"
echo \
"chr1	10	20	a1
chr1	30	100	a2,a3,a4" > exp
$BT merge -i a.names.bed -c 4 -o collapse > obs
check obs exp
rm obs exp

###########################################################
# collapsed list of the names, and sum of the scores
###########################################################
echo "    merge.t6...\c"
echo \
"chr1	10	20	a1	1
chr1	30	100	a2,a3,a4	9
chr2	10	20	a1	5
chr2	30	40	a2	6
chr2	42	100	a3,a4	15" > exp
$BT merge -i a.full.bed -c 4,5  -o collapse,sum > obs
check obs exp
rm obs exp

###########################################################
# count intervals and sum of scores
###########################################################
echo "    merge.t7...\c"
echo \
"chr1	10	20	1	1
chr1	30	100	3	9
chr2	10	20	1	5
chr2	30	40	1	6
chr2	42	100	2	15" > exp
$BT merge -i a.full.bed -c 5 -o count,sum> obs
check obs exp
rm obs exp

###########################################################
# count, collapsed names, and sum of scores
###########################################################
echo "    merge.t8...\c"
echo \
"chr1	10	20	a1	1	1
chr1	30	100	a2,a3,a4	9	3
chr2	10	20	a1	5	1
chr2	30	40	a2	6	1
chr2	42	100	a3,a4	15	2" > exp
$BT merge -i a.full.bed -c 4,5,4 -o collapse,sum,count > obs
check obs exp
rm obs exp

###########################################################
# stranded merge, show sign, collapsed names, sum of
# scores, and count
###########################################################
echo "    merge.t9...\c"
echo \
"chr1	10	20	+	a1	1	1
chr1	30	40	+	a2	2	1
chr1	40	50	-	a3	3	1
chr1	45	100	+	a4	4	1
chr2	10	20	+	a1	5	1
chr2	30	40	+	a2	6	1
chr2	42	50	+	a3	7	1
chr2	45	100	-	a4	8	1" > exp
$BT merge -i a.full.bed -s -c 6,4,5,6 -o distinct,collapse,sum,count > obs
check obs exp
rm obs exp

###########################################################
#  Test the use of a custom delimiter for -nms
###########################################################
echo "    merge.t10...\c"
echo \
"chr1	10	20	a1
chr1	30	100	a2|a3|a4" > exp
$BT merge -i a.names.bed -delim "|" -c 4 -o collapse > obs
check obs exp
rm obs exp
