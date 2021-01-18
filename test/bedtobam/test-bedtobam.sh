set -e;
BT=${BT-../../bin/bedtools}
htsutil=${htsutil-../htsutil}

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

##################################################################
#  Test one block without -split
##################################################################
echo -e "    bedtobam.t1...\c"
echo \
"read_name	0	1	1001	255	1000M	*	0	0	*	*">exp
echo -e "1\t1000\t2000\tread_name\t255\t+" | $BT bedtobam -i - -g chrsize.tmp| $htsutil viewbamrecords > obs
check obs exp
rm obs exp

[[ $FAILURES -eq 0 ]] || exit 1;
