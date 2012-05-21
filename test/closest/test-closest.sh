BT=../../bin/bedtools

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
