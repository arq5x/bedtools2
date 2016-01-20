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
# test shifting forward via -s
###########################################################
echo "    shift.t1...\c"
echo \
"chr1	105	205	a1	1	+
chr1	105	205	a2	2	-" > exp
$BT shift -i a.bed -s 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test shifting backward via -s
###########################################################
echo "    shift.t2...\c"
echo \
"chr1	95	195	a1	1	+
chr1	95	195	a2	2	-" > exp
$BT shift -i a.bed -s -5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test shifting forward via -m and -p
###########################################################
echo "    shift.t3...\c"
echo \
"chr1	105	205	a1	1	+
chr1	105	205	a2	2	-" > exp
$BT shift -i a.bed -p 5 -m 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test shifting backward via -m and -p
###########################################################
echo "    shift.t3...\c"
echo \
"chr1	95	195	a1	1	+
chr1	95	195	a2	2	-" > exp
$BT shift -i a.bed -p -5 -m -5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -m shift (-p == 0)
###########################################################
echo "    shift.t4...\c"
echo \
"chr1	100	200	a1	1	+
chr1	95	195	a2	2	-" > exp
$BT shift -i a.bed -m -5 -p 0 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test just a -p shift (-m == 0)
###########################################################
echo "    shift.t5...\c"
echo \
"chr1	105	205	a1	1	+
chr1	100	200	a2	2	-" > exp
$BT shift -i a.bed -m 0 -p 5 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the start of the chrom
###########################################################
echo "    shift.t6...\c"
echo \
"chr1	0	1	a1	1	+
chr1	0	1	a2	2	-" > exp
$BT shift -i a.bed -s -200 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test going beyond the end of the chrom
###########################################################
echo "    shift.t7...\c"
echo \
"chr1	999	1000	a1	1	+
chr1	999	1000	a2	2	-" > exp
$BT shift -i a.bed -s 1000 -g tiny.genome > obs
check obs exp
rm obs exp

###########################################################
# test shift being larger than a signed int
###########################################################
echo "    shift.t8...\c"
echo \
"chr1	999	1000	a1	1	+
chr1	999	1000	a2	2	-" > exp
$BT shift -i a.bed -s 3000000000 -g tiny.genome > obs
check obs exp
rm obs exp
