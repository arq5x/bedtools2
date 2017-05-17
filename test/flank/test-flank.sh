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
# chr1	100	200	a1	1	+
# chr1	100	200	a2	2	-

###########################################################
# test matching flanks via -b
###########################################################
echo -e "    flank.t1...\c"
echo \
"chr1	95	100	a1	1	+
chr1	200	205	a1	1	+
chr1	95	100	a2	2	-
chr1	200	205	a2	2	-" > exp
$BT flank -i a.bed -b 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test matching flanks via -l and -r
###########################################################
echo -e "    flank.t2...\c"
echo \
"chr1	95	100	a1	1	+
chr1	200	205	a1	1	+
chr1	95	100	a2	2	-
chr1	200	205	a2	2	-" > exp
$BT flank -i a.bed -l 5 -r 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -l flank (-r == 0)
###########################################################
echo -e "    flank.t3...\c"
echo \
"chr1	95	100	a1	1	+
chr1	95	100	a2	2	-" > exp
$BT flank -i a.bed -l 5 -r 0 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -r flank (-l == 0)
###########################################################
echo -e "    flank.t4...\c"
echo \
"chr1	200	205	a1	1	+
chr1	200	205	a2	2	-" > exp
$BT flank -i a.bed -l 0 -r 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -l flank (-r == 0) with -s
###########################################################
echo -e "    flank.t5...\c"
echo \
"chr1	95	100	a1	1	+
chr1	200	205	a2	2	-" > exp
$BT flank -i a.bed -l 5 -r 0 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -r flank (-l == 0) with -s
###########################################################
echo -e "    flank.t6...\c"
echo \
"chr1	200	205	a1	1	+
chr1	95	100	a2	2	-" > exp
$BT flank -i a.bed -l 0 -r 5 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test -b with -s
###########################################################
echo -e "    flank.t7...\c"
echo \
"chr1	95	100	a1	1	+
chr1	200	205	a1	1	+
chr1	95	100	a2	2	-
chr1	200	205	a2	2	-" > exp
$BT flank -i a.bed -b 5 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start of the chrom
###########################################################
echo -e "    flank.t8...\c"
echo \
"chr1	0	100	a1	1	+
chr1	200	400	a1	1	+
chr1	0	100	a2	2	-
chr1	200	400	a2	2	-" > exp
$BT flank -i a.bed -b 200 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the end of the chrom
###########################################################
echo -e "    flank.t9...\c"
echo \
"chr1	200	1000	a1	1	+
chr1	200	1000	a2	2	-" > exp
$BT flank -i a.bed -l 0 -r 1000 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start and end of the chrom
###########################################################
echo -e "    flank.t10...\c"
echo \
"chr1	0	100	a1	1	+
chr1	200	1000	a1	1	+
chr1	0	100	a2	2	-
chr1	200	1000	a2	2	-" > exp
$BT flank -i a.bed -b 2000 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start and end of the chrom with -s
###########################################################
echo -e "    flank.t11...\c"
echo \
"chr1	0	100	a1	1	+
chr1	200	1000	a1	1	+
chr1	0	100	a2	2	-
chr1	200	1000	a2	2	-" > exp
$BT flank -i a.bed -b 2000 -s -g tiny.genome > obs
check obs exp
rm obs exp
[[ $FAILURES -eq 0 ]] || exit 1;
