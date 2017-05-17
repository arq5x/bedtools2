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

# cat in.bed
#chr1	10	20	1,2,3	10,20,30
#chr1	40	50	4,5,6	40,50,60

###########################################################
#  expand one column
###########################################################
echo -e "    expand.t1...\c"
echo \
"chr1	10	20	1	10,20,30
chr1	10	20	2	10,20,30
chr1	10	20	3	10,20,30
chr1	40	50	4	40,50,60
chr1	40	50	5	40,50,60
chr1	40	50	6	40,50,60" > exp
$BT expand -i expand.txt -c 4 > obs
check obs exp
rm obs exp

###########################################################
#  expand multiple columns
###########################################################
echo -e "    expand.t2...\c"
echo \
"chr1	10	20	1	10
chr1	10	20	2	20
chr1	10	20	3	30
chr1	40	50	4	40
chr1	40	50	5	50
chr1	40	50	6	60" > exp
$BT expand -i expand.txt -c 4,5 > obs
check obs exp
rm obs exp

###########################################################
#  expand multiple columns while switching order
###########################################################
echo -e "    expand.t3...\c"
echo \
"chr1	10	20	10	1
chr1	10	20	20	2
chr1	10	20	30	3
chr1	40	50	40	4
chr1	40	50	50	5
chr1	40	50	60	6" > exp
$BT expand -i expand.txt -c 5,4 > obs
check obs exp
rm obs exp
[[ $FAILURES -eq 0 ]] || exit 1;
