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
# chr1	100	200	a1	1	+
# chr1	100	200	a2	2	-

###########################################################
# test matching flanks via -b
###########################################################
echo "    slop.t1...\c"
echo \
"chr1	95	205	a1	1	+
chr1	95	205	a2	2	-" > exp
$BT slop -i a.bed -b 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test matching flanks via -l and -r
###########################################################
echo "    slop.t2...\c"
echo \
"chr1	95	205	a1	1	+
chr1	95	205	a2	2	-" > exp > exp
$BT slop -i a.bed -l 5 -r 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -l flank (-r == 0)
###########################################################
echo "    slop.t3...\c"
echo \
"chr1	95	200	a1	1	+
chr1	95	200	a2	2	-" > exp
$BT slop -i a.bed -l 5 -r 0 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -r flank (-l == 0)
###########################################################
echo "    slop.t4...\c"
echo \
"chr1	100	205	a1	1	+
chr1	100	205	a2	2	-" > exp
$BT slop -i a.bed -l 0 -r 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -l flank (-r == 0) with -s
###########################################################
echo "    slop.t5...\c"
echo \
"chr1	95	200	a1	1	+
chr1	100	205	a2	2	-" > exp
$BT slop -i a.bed -l 5 -r 0 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -r flank (-l == 0) with -s
###########################################################
echo "    slop.t6...\c"
echo \
"chr1	100	205	a1	1	+
chr1	95	200	a2	2	-" > exp
$BT slop -i a.bed -l 0 -r 5 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test -b with -s
###########################################################
echo "    slop.t7...\c"
echo \
"chr1	95	205	a1	1	+
chr1	95	205	a2	2	-" > exp
$BT slop -i a.bed -b 5 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start of the chrom
###########################################################
echo "    slop.t8...\c"
echo \
"chr1	0	400	a1	1	+
chr1	0	400	a2	2	-" > exp
$BT slop -i a.bed -b 200 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the end of the chrom
###########################################################
echo "    slop.t9...\c"
echo \
"chr1	100	1000	a1	1	+
chr1	100	1000	a2	2	-" > exp
$BT slop -i a.bed -l 0 -r 1000 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start and end of the chrom
###########################################################
echo "    slop.t10...\c"
echo \
"chr1	0	1000	a1	1	+
chr1	0	1000	a2	2	-" > exp
$BT slop -i a.bed -b 2000 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start and end of the chrom with -s
###########################################################
echo "    slop.t11...\c"
echo \
"chr1	0	1000	a1	1	+
chr1	0	1000	a2	2	-" > exp
$BT slop -i a.bed -b 2000 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test slop factor being larger than a signed int
###########################################################
echo "    slop.t12...\c"
echo \
"chr1	0	1000	a1	1	+
chr1	0	1000	a2	2	-" > exp
$BT slop -i a.bed -b 3000000000 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test that old floating-point issues are solved
###########################################################
echo "    slop.t13...\c"
echo \
"chr1	16778071	16778771" > exp
echo -e "chr1\t16778271\t16778571"| $BT slop -l 200 -r 200 -i - -g ../../genomes/human.hg19.genome > obs
check obs exp
rm obs exp

###########################################################
# test that old floating-point issues are solved
###########################################################
echo "    slop.t14...\c"
echo \
"chr1	16778072	16778772" > exp
echo -e "chr1\t16778272\t16778572"| $BT slop -l 200 -r 200 -i - -g ../../genomes/human.hg19.genome > obs
check obs exp
rm obs exp

echo "    slop.t15...\c"
echo \
"chr1	159	171
chr1	90	210" > exp
echo -e "chr1\t160\t170\nchr1\t100\t200"| $BT slop -b 0.1 -pct -i - -g ../../genomes/human.hg19.genome > obs
check obs exp
rm obs exp
