set -e;
BT="${BT-../../bin/bedtools}"

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
# chr1	100	200	a1	1	+
# chr1	100	200	a2	2	-

###########################################################
# test simple
###########################################################
echo -e "    revcomp.t1...\c"
echo \
"chr1	800	900	a1	1	+
chr1	800	900	a2	2	-" > exp
"$BT" revcomp -i a.bed -g tiny.genome > obs
check obs exp
rm obs exp

echo -e "    revcomp.t2...\c"
echo \
"chr1	249250421	249250521	a1	1	+
chr1	249250421	249250521	a2	2	-" > exp
"$BT" revcomp -i a.bed -g huge.genome > obs
check obs exp
rm obs exp


###########################################################
# test clamp
###########################################################
echo -e "    revcomp.t3...\c"
echo \
"chr1	0	0	NM_032291	0	+
chr1	0	0	NR_036634	0	-" > exp
"$BT" revcomp -i b.bed -g tiny.genome > obs
check obs exp
rm obs exp


[[ $FAILURES -eq 0 ]] || exit 1;
