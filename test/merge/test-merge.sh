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
# Test #1
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
#
# NOTE: Testing for sorted input is now deprecated, as the
# FileRecordMgr is already testing for that.
#
###########################################################
# Test #2
#  Enforce coordinate sorted input.
###########################################################
#echo "    merge.t2...\c"
#command -v tac 2>/dev/null || alias tac="sed '1!G;h;\$!d'"
#tac a.bed | $BT merge -i - 2> obs
#echo "ERROR: input file: (-) is not sorted by chrom then start.
#       The start coordinate at line 3 is less than the start at line 2" > exp
#check obs exp
#rm obs exp


###########################################################
# Test #3
#  Test the counting of merged intervals. (-n)
###########################################################
echo "    merge.t3...\c"
echo \
"chr1	10	20	1
chr1	30	100	3" > exp
$BT merge -i a.bed -n > obs
check obs exp
rm obs exp


###########################################################
# Test #4
#  Test the listing of names from merged intervals. (-nms)
#  a.bed should fail, as there is no name field
###########################################################
echo "    merge.t4...\c"
echo \
"*****
***** ERROR: Requested column 4, but database file a.bed only has fields 1 - 3." > exp
$BT merge -i a.bed -nms 2>&1 > /dev/null | head -3 | tail -2 > obs
check obs exp
rm obs exp


###########################################################
# Test #5
#  Test the listing of names from merged intervals. (-nms)
#  a.named.bed should work, as there are name fields
#  
# cat a.names.bed
# chr1	10	20	a1
# chr1	30	40	a2
# chr1	40	50	a3
# chr1	45	100	a4
###########################################################
echo "    merge.t5...\c"
echo \
"chr1	10	20	a1
chr1	30	100	a2,a3,a4" > exp
$BT merge -i a.names.bed -nms > obs
check obs exp
rm obs exp

###########################################################
# -nms and -scores sum
###########################################################
echo "    merge.t6...\c"
echo \
"chr1	10	20	a1	1
chr1	30	100	a2,a3,a4	9
chr2	10	20	a1	5
chr2	30	40	a2	6
chr2	42	100	a3,a4	15" > exp
$BT merge -i a.full.bed -nms -scores sum> obs
check obs exp
rm obs exp

###########################################################
# -n and -scores sum
###########################################################
echo "    merge.t7...\c"
echo \
"chr1	10	20	1	1
chr1	30	100	3	9
chr2	10	20	1	5
chr2	30	40	1	6
chr2	42	100	2	15" > exp
$BT merge -i a.full.bed -n -scores sum> obs
check obs exp
rm obs exp

###########################################################
# -n, -nms, and -scores sum
###########################################################
echo "    merge.t8...\c"
echo \
"chr1	10	20	a1	1	1
chr1	30	100	a2,a3,a4	9	3
chr2	10	20	a1	5	1
chr2	30	40	a2	6	1
chr2	42	100	a3,a4	15	2" > exp
$BT merge -i a.full.bed -nms -scores sum -n> obs
check obs exp
rm obs exp

###########################################################
# -s, -n, -nms, and -scores sum
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
$BT merge -i a.full.bed -s -nms -scores sum -n> obs
check obs exp
rm obs exp

###########################################################
# Test #10
#  Test the use of a custom delimiter for -nms
#  
# cat a.names.bed
# chr1	10	20	a1
# chr1	30	40	a2
# chr1	40	50	a3
# chr1	45	100	a4
###########################################################
echo "    merge.t10...\c"
echo \
"chr1	10	20	a1
chr1	30	100	a2|a3|a4" > exp
$BT merge -i a.names.bed -nms -delim "|" > obs
check obs exp
rm obs exp
