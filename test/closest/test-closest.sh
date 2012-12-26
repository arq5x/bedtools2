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
# chr1	20	21

# cat b.bed
# chr1	20	21

# cat b-one-bp-closer.bed
# chr1	19	21

###########################################################
# test 1bp apart; checking for off-by-one errors
###########################################################
echo "    closest.t1...\c"
echo \
"chr1	10	20	chr1	20	21	1" > exp
$BT closest -a a.bed -b b.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test reciprocal of t1
###########################################################
echo "    closest.t2...\c"
echo \
"chr1	20	21	chr1	10	20	1" > exp
$BT closest -a b.bed -b a.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test 0bp apart; checking for off-by-one errors
###########################################################
echo "    closest.t3...\c"
echo \
"chr1	10	20	chr1	19	21	0" > exp
$BT closest -a a.bed -b b-one-bp-closer.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test reciprocal of t3
###########################################################
echo "    closest.t4...\c"
echo \
"chr1	19	21	chr1	10	20	0" > exp
$BT closest -a b-one-bp-closer.bed -b a.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test closest without forcing different names ( -N )
###########################################################
echo "    closest.t5...\c"
echo \
"chr1	10	20	break1	chr1	40	50	break1	21
chr1	55	58	break2	chr1	60	70	break2	3" > exp
$BT closest -a a.names.bed -b b.names.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test closest with forcing different names ( -N )
###########################################################
echo "    closest.t6...\c"
echo \
"chr1	10	20	break1	chr1	60	70	break2	41
chr1	55	58	break2	chr1	40	50	break1	6" > exp
$BT closest -a a.names.bed -b b.names.bed -d -N > obs
check obs exp
rm obs exp

###########################################################
# test closest forcing -s yet no matching strands on chrom
###########################################################
echo "    closest.t7...\c"
echo \
"chr1	100	200	a	10	+	.	-1	-1	.	-1	." > exp
$BT closest -a strand-test-a.bed -b strand-test-b.bed -s > obs
check obs exp
rm obs exp

###########################################################
# test closest forcing -S with only an opp strands on chrom
###########################################################
echo "    closest.t8...\c"
echo \
"chr1	100	200	a	10	+	chr1	90	120	b	1	-" > exp
$BT closest -a strand-test-a.bed -b strand-test-b.bed -S > obs
check obs exp
rm obs exp


###########################################################
# test reproting of all overlapping features
###########################################################
echo "    closest.t9...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	150	201
chr1	200	201	chr1	175	375
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed > obs
check obs exp
rm obs exp

###########################################################
# test reproting of first overlapping feature
###########################################################
echo "    closest.t10...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	150	201
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed -t first > obs
check obs exp
rm obs exp

###########################################################
# test reproting of last overlapping feature
###########################################################
echo "    closest.t11...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	175	375
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed -t last > obs
check obs exp
rm obs exp

###########################################################
# test reproting of no overlapping feature
###########################################################
echo "    closest.t12...\c"
echo \
"chr1	100	101	chr1	150	201
chr1	200	201	chr1	100	101
chr1	300	301	chr1	150	201
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed -io > obs
check obs exp
rm obs exp
