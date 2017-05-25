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
echo -e "    slop.t1...\c"
echo \
"chr1	95	205	a1	1	+
chr1	95	205	a2	2	-" > exp
$BT slop -i a.bed -b 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test matching flanks via -l and -r
###########################################################
echo -e "    slop.t2...\c"
echo \
"chr1	95	205	a1	1	+
chr1	95	205	a2	2	-" > exp > exp
$BT slop -i a.bed -l 5 -r 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -l flank (-r == 0)
###########################################################
echo -e "    slop.t3...\c"
echo \
"chr1	95	200	a1	1	+
chr1	95	200	a2	2	-" > exp
$BT slop -i a.bed -l 5 -r 0 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -r flank (-l == 0)
###########################################################
echo -e "    slop.t4...\c"
echo \
"chr1	100	205	a1	1	+
chr1	100	205	a2	2	-" > exp
$BT slop -i a.bed -l 0 -r 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -l flank (-r == 0) with -s
###########################################################
echo -e "    slop.t5...\c"
echo \
"chr1	95	200	a1	1	+
chr1	100	205	a2	2	-" > exp
$BT slop -i a.bed -l 5 -r 0 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -r flank (-l == 0) with -s
###########################################################
echo -e "    slop.t6...\c"
echo \
"chr1	100	205	a1	1	+
chr1	95	200	a2	2	-" > exp
$BT slop -i a.bed -l 0 -r 5 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test -b with -s
###########################################################
echo -e "    slop.t7...\c"
echo \
"chr1	95	205	a1	1	+
chr1	95	205	a2	2	-" > exp
$BT slop -i a.bed -b 5 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start of the chrom
###########################################################
echo -e "    slop.t8...\c"
echo \
"chr1	0	400	a1	1	+
chr1	0	400	a2	2	-" > exp
$BT slop -i a.bed -b 200 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the end of the chrom
###########################################################
echo -e "    slop.t9...\c"
echo \
"chr1	100	1000	a1	1	+
chr1	100	1000	a2	2	-" > exp
$BT slop -i a.bed -l 0 -r 1000 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start and end of the chrom
###########################################################
echo -e "    slop.t10...\c"
echo \
"chr1	0	1000	a1	1	+
chr1	0	1000	a2	2	-" > exp
$BT slop -i a.bed -b 2000 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start and end of the chrom with -s
###########################################################
echo -e "    slop.t11...\c"
echo \
"chr1	0	1000	a1	1	+
chr1	0	1000	a2	2	-" > exp
$BT slop -i a.bed -b 2000 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test slop factor being larger than a signed int
###########################################################
echo -e "    slop.t12...\c"
echo \
"chr1	0	1000	a1	1	+
chr1	0	1000	a2	2	-" > exp
$BT slop -i a.bed -b 3000000000 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test that old floating-point issues are solved
###########################################################
echo -e "    slop.t13...\c"
echo \
"chr1	16778071	16778771" > exp
echo -e "chr1\t16778271\t16778571"| $BT slop -l 200 -r 200 -i - -g ../../genomes/human.hg19.genome > obs
check obs exp
rm obs exp

###########################################################
# test that old floating-point issues are solved
###########################################################
echo -e "    slop.t14...\c"
echo \
"chr1	16778072	16778772" > exp
echo -e "chr1\t16778272\t16778572"| $BT slop -l 200 -r 200 -i - -g ../../genomes/human.hg19.genome > obs
check obs exp
rm obs exp

echo -e "    slop.t15...\c"
echo \
"chr1	159	171
chr1	90	210" > exp
echo -e "chr1\t160\t170\nchr1\t100\t200"| $BT slop -b 0.1 -pct -i - -g ../../genomes/human.hg19.genome > obs
check obs exp
rm obs exp
###########################################################
# test negative slop on l with strand
###########################################################

echo -e "    slop.t16...\c"
echo \
"chr1	360	380" > exp
echo -e "chr1\t300\t320" | $BT slop -l -60 -r 60 -i - -g tiny.genome > obs
check obs exp
rm obs exp

echo -e "    slop.t17...\c"
echo \
"chr1	240	260	a1	5	-" > exp
echo -e "chr1\t300\t320\ta1\t5\t-" | $BT slop -s -l -60 -r 60 -i - -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test negative slop on r with strand
###########################################################

echo -e "    slop.t18...\c"
echo \
"chr1	240	260" > exp
echo -e "chr1\t300\t320" | $BT slop -l 60 -r -60 -i - -g tiny.genome > obs
check obs exp
rm obs exp

echo -e "    slop.t19...\c"
echo \
"chr1	360	380	a1	5	-" > exp
echo -e "chr1\t300\t320\ta1\t5\t-" | $BT slop -s -l 60 -r -60 -i - -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test crossover during negative slop
###########################################################

echo -e "    slop.t20...\c"
echo \
"chr1	260	360	a1	5	-" > exp
echo -e "chr1\t300\t320\ta1\t5\t-" | $BT slop -s -l -60 -r -60 -i - -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test edge cases
###########################################################

echo -e "    slop.t21...\c"
echo \
"chr1	999	1000	a1	5	-" > exp
echo -e "chr1\t950\t970\ta1\t5\t-" | $BT slop -s -l 60 -r -60 -i - -g tiny.genome > obs
check obs exp
rm obs exp

echo -e "    slop.t22...\c"
echo \
"chr1	0	1	a1	5	-" > exp
echo -e "chr1\t50\t60\ta1\t5\t-" | $BT slop -s -l -60 -r 60 -i - -g tiny.genome > obs
check obs exp
rm obs exp

[[ $FAILURES -eq 0 ]] || exit 1;
