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
# chr1	100	200	a1	1	+
# chr1	100	200	a2	2	-

###########################################################
# test matching flanks via -b
###########################################################
echo "    flank.t1...\c"
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
echo "    flank.t2...\c"
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
echo "    flank.t3...\c"
echo \
"chr1	95	100	a1	1	+
chr1	95	100	a2	2	-" > exp
$BT flank -i a.bed -l 5 -r 0 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -r flank (-l == 0)
###########################################################
echo "    flank.t4...\c"
echo \
"chr1	200	205	a1	1	+
chr1	200	205	a2	2	-" > exp
$BT flank -i a.bed -l 0 -r 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -l flank (-r == 0) with -s
###########################################################
echo "    flank.t5...\c"
echo \
"chr1	95	100	a1	1	+
chr1	200	205	a2	2	-" > exp
$BT flank -i a.bed -l 5 -r 0 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -r flank (-l == 0) with -s
###########################################################
echo "    flank.t6...\c"
echo \
"chr1	200	205	a1	1	+
chr1	95	100	a2	2	-" > exp
$BT flank -i a.bed -l 0 -r 5 -s -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test -b with -s
###########################################################
echo "    flank.t7...\c"
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
echo "    flank.t8...\c"
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
echo "    flank.t9...\c"
echo \
"chr1	200	1000	a1	1	+
chr1	200	1000	a2	2	-" > exp
$BT flank -i a.bed -l 0 -r 1000 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start and end of the chrom
###########################################################
echo "    flank.t10...\c"
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
echo "    flank.t11...\c"
echo \
"chr1	0	100	a1	1	+
chr1	200	1000	a1	1	+
chr1	0	100	a2	2	-
chr1	200	1000	a2	2	-" > exp
$BT flank -i a.bed -b 2000 -s -g tiny.genome > obs
check obs exp
rm obs exp