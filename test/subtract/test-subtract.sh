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
# chr1	10	20	a1	1	+
# chr1	50	70	a2	2	-
# 
# cat b.bed
# chr1	18	25	b1	1	-
# chr1	80	90	b2	2	+


###########################################################
# test baseline subtraction
###########################################################
echo -e "    subtract.t1...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo -e "    subtract.t2...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.5 > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo -e "    subtract.t3...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.1 > obs
check obs exp
rm obs exp


###########################################################
# test -s subtraction
###########################################################
echo -e "    subtract.t4...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -s > obs
check obs exp
rm obs exp


###########################################################
# test -S subtraction
###########################################################
echo -e "    subtract.t5...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -S > obs
check obs exp
rm obs exp


###########################################################
# test -A subtraction
###########################################################
echo -e "    subtract.t6...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A > obs
check obs exp
rm obs exp


###########################################################
# test -A with -f subtraction
###########################################################
echo -e "    subtract.t7...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.5 > obs
check obs exp
rm obs exp

###########################################################
# test -A with -f subtraction
###########################################################
echo -e "    subtract.t8...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.1 > obs
check obs exp
rm obs exp


###########################################################
# test -N with -f subtraction
###########################################################
echo -e "    subtract.t9...\c"
echo \
"chr1	0	10" > exp
$BT subtract -a c.bed -b d.bed -N -f 0.4 > obs
check obs exp
rm obs exp

###########################################################
# test -N with -f subtraction
###########################################################
echo -e "    subtract.t10...\c"
touch exp
$BT subtract -a c.bed -b d.bed -N -f 0.39 > obs
check obs exp
rm obs exp


###########################################################
# 
# REPEAT TESTS 1-8 WITH -sorted 
#
###########################################################


###########################################################
# test baseline subtraction
###########################################################
echo -e "    subtract.t11...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo -e "    subtract.t12...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.5  -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -f subtraction
###########################################################
echo -e "    subtract.t13...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -f 0.1  -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -s subtraction
###########################################################
echo -e "    subtract.t14...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -s -sorted  > obs
check obs exp
rm obs exp


###########################################################
# test -S subtraction
###########################################################
echo -e "    subtract.t15...\c"
echo \
"chr1	10	18	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -S -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -A subtraction
###########################################################
echo -e "    subtract.t16...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A  -sorted > obs
check obs exp
rm obs exp


###########################################################
# test -A with -f subtraction
###########################################################
echo -e "    subtract.t17...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.5  -sorted > obs
check obs exp
rm obs exp

###########################################################
# test -A with -f subtraction
###########################################################
echo -e "    subtract.t18...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed -A -f 0.1  -sorted > obs
check obs exp
rm obs exp


###########################################################
#
# TEST WITH MULTIPLE DATABASES
# 
###########################################################


###########################################################
# test with 2 DBs
###########################################################
echo -e "    subtract.t19...\c"
echo \
"chr1	15	18	a1	1	+
chr1	50	55	a2	2	-
chr1	65	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed b2.bed > obs
check obs exp
rm obs exp

###########################################################
# test with 2 DBs, -f option
###########################################################
echo -e "    subtract.t20...\c"
echo \
"chr1	10	20	a1	1	+
chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed b2.bed -f 0.8 > obs
check obs exp
rm obs exp

###########################################################
# test with 2 DBs, -f option and -N options
###########################################################
echo -e "    subtract.t21...\c"
echo \
"chr1	50	70	a2	2	-" > exp
$BT subtract -a a.bed -b b.bed b2.bed -f 0.6 -N > obs
check obs exp
rm obs exp

###########################################################
###########################################################
#  subtract.new.*: tests added for the intersect rewrite.
#  subtract runs on intersect's readers, bin index and
#  overlap test (ContextSubtract subclasses ContextIntersect),
#  so each test pins a visible output of that machinery.
###########################################################
###########################################################
htsutil=${htsutil-../htsutil}

###########################################################
#  subtract.new.t001: B over the left end of A trims A's start
#  A chr1:100-200; B 50-120 removes bases 100-119, leaving
#  120-200 with A's other fields unchanged.
###########################################################
echo -e "    subtract.new.t001...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t50\t120\n' > b.tmp.bed
printf 'chr1\t120\t200\ta1\t0\t+\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t002: B over the right end of A trims A's end
#  B 180-250 removes bases 180-199, leaving 100-180.
###########################################################
echo -e "    subtract.new.t002...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t180\t250\n' > b.tmp.bed
printf 'chr1\t100\t180\ta1\t0\t+\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t003: B inside A splits A into two pieces
#  B 140-160 leaves 100-140 and 160-200, each a copy of A's
#  other fields.
###########################################################
echo -e "    subtract.new.t003...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t140\t160\n' > b.tmp.bed
printf 'chr1\t100\t140\ta1\t0\t+\nchr1\t160\t200\ta1\t0\t+\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t004: B identical to A removes A entirely
###########################################################
echo -e "    subtract.new.t004...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t005: B containing A removes A entirely
###########################################################
echo -e "    subtract.new.t005...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t50\t300\n' > b.tmp.bed
printf '' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t006: three disjoint B split A into four pieces
#  B is listed out of order (170, 110, 140); the pieces are
#  100-110, 120-140, 150-170, 180-200 in ascending order.
###########################################################
echo -e "    subtract.new.t006...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
printf 'chr1\t170\t180\nchr1\t110\t120\nchr1\t140\t150\n' > b.tmp.bed
printf 'chr1\t100\t110\ta1\nchr1\t120\t140\ta1\nchr1\t150\t170\ta1\nchr1\t180\t200\ta1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t007: overlapping, nested and duplicated B
#  B 110-130, 120-140, 125-128 (nested), 110-130 (duplicate),
#  160-170, 165-175. Their union is 110-140 and 160-175, so
#  100-110, 140-160 and 175-200 remain.
###########################################################
echo -e "    subtract.new.t007...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t110\t130\nchr1\t120\t140\nchr1\t125\t128\nchr1\t110\t130\nchr1\t160\t170\nchr1\t165\t175\n' > b.tmp.bed
printf 'chr1\t100\t110\nchr1\t140\t160\nchr1\t175\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t008: B touching either end of A removes nothing
#  BED is half-open: B 200-300 starts where A ends and B 50-100
#  ends where A starts, so neither shares a base with A.
###########################################################
echo -e "    subtract.new.t008...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t200\t300\nchr1\t50\t100\n' > b.tmp.bed
printf 'chr1\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t009: 1 bp B on A's first and last base
#  B 100-101 and 199-200 remove one base at each end: 101-199.
###########################################################
echo -e "    subtract.new.t009...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t101\nchr1\t199\t200\n' > b.tmp.bed
printf 'chr1\t101\t199\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t010: 1 bp A records
#  x1 (100-101) is covered by B 100-101 and removed; x2
#  (200-201) only touches B 201-202 and is kept.
###########################################################
echo -e "    subtract.new.t010...\c"
printf 'chr1\t100\t101\tx1\nchr1\t200\t201\tx2\n' > a.tmp.bed
printf 'chr1\t100\t101\nchr1\t201\t202\n' > b.tmp.bed
printf 'chr1\t200\t201\tx2\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t011: pieces 1 bp long
#  A 100-105 minus 101-102 and 103-104 leaves 100-101,
#  102-103 and 104-105.
###########################################################
echo -e "    subtract.new.t011...\c"
printf 'chr1\t100\t105\n' > a.tmp.bed
printf 'chr1\t101\t102\nchr1\t103\t104\n' > b.tmp.bed
printf 'chr1\t100\t101\nchr1\t102\t103\nchr1\t104\t105\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t012: B at every second base gives 50 pieces
#  A 100-200; B holds the 50 one-base intervals at even
#  positions 100, 102, ..., 198, all in one bin. The odd bases
#  remain, as 50 one-base pieces in ascending order.
###########################################################
echo -e "    subtract.new.t012...\c"
printf 'chr1\t100\t200\tw\n' > a.tmp.bed
i=100
: > b.tmp.bed
: > exp
while [ $i -lt 200 ]; do
	printf 'chr1\t%d\t%d\n' $i $((i + 1)) >> b.tmp.bed
	printf 'chr1\t%d\t%d\tw\n' $((i + 1)) $((i + 2)) >> exp
	i=$((i + 2))
done
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t013: duplicated A records are each reported
#  Both copies of A lose 150-160, and the pieces follow A's
#  input order: first copy's two pieces, then the second's.
###########################################################
echo -e "    subtract.new.t013...\c"
printf 'chr1\t100\t200\td\nchr1\t100\t200\td\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t150\td\nchr1\t160\t200\td\nchr1\t100\t150\td\nchr1\t160\t200\td\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t014: overlapping A records are not merged
#  r1 100-200 and r2 150-250 each lose B 180-190 on their own:
#  r1 -> 100-180, 190-200; r2 -> 150-180, 190-250.
###########################################################
echo -e "    subtract.new.t014...\c"
printf 'chr1\t100\t200\tr1\nchr1\t150\t250\tr2\n' > a.tmp.bed
printf 'chr1\t180\t190\n' > b.tmp.bed
printf 'chr1\t100\t180\tr1\nchr1\t190\t200\tr1\nchr1\t150\t180\tr2\nchr1\t190\t250\tr2\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t015: unsorted A is reported in input order
#  A is out of order by chromosome and position. Without
#  -sorted each A record is processed as read: q1 and q2 are
#  untouched, q3 loses 150-160, q4 loses 5-10.
###########################################################
echo -e "    subtract.new.t015...\c"
printf 'chr2\t50\t60\tq1\nchr1\t300\t400\tq2\nchr1\t100\t200\tq3\nchr10\t5\t15\tq4\n' > a.tmp.bed
printf 'chr1\t150\t160\nchr10\t0\t10\n' > b.tmp.bed
printf 'chr2\t50\t60\tq1\nchr1\t300\t400\tq2\nchr1\t100\t150\tq3\nchr1\t160\t200\tq3\nchr10\t10\t15\tq4\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t016: chromosomes present in only one file
#  chr1 is only in A, so its record is kept whole; chr2 and
#  chrY are only in B and add nothing; chrX loses 150-160.
###########################################################
echo -e "    subtract.new.t016...\c"
printf 'chr1\t100\t200\nchrX\t100\t200\n' > a.tmp.bed
printf 'chr2\t100\t200\nchrX\t150\t160\nchrY\t0\t1000\n' > b.tmp.bed
printf 'chr1\t100\t200\nchrX\t100\t150\nchrX\t160\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t017: chromosome names must match exactly
#  B on chr10, chr11 and chr1_alt covers the same coordinates
#  as A on chr1 but is on other chromosomes, so A is kept.
###########################################################
echo -e "    subtract.new.t017...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr10\t100\t200\nchr11\t100\t200\nchr1_alt\t100\t200\n' > b.tmp.bed
printf 'chr1\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t018: records straddling bin boundaries
#  The finest bins are 16384 bp, the next 131072 bp, then
#  1048576 bp. e1-e3 and their B straddle those edges; e4 lies
#  in one fine bin while its B straddles 16384; e5 ends at
#  131072 while its B straddles it.
#  e1 16380-16390 minus 16383-16385 -> 16380-16383, 16385-16390
#  e2 131070-131080 minus 131071-131073 -> 131070-131071, 131073-131080
#  e3 1048570-1048580 minus 1048575-1048577 -> 1048570-1048575, 1048577-1048580
#  e4 16384-16390 minus 16383-16385 -> 16385-16390
#  e5 131060-131072 minus 131071-131073 -> 131060-131071
###########################################################
echo -e "    subtract.new.t018...\c"
printf 'chr1\t16380\t16390\te1\nchr1\t131070\t131080\te2\nchr1\t1048570\t1048580\te3\nchr1\t16384\t16390\te4\nchr1\t131060\t131072\te5\n' > a.tmp.bed
printf 'chr1\t16383\t16385\nchr1\t131071\t131073\nchr1\t1048575\t1048577\n' > b.tmp.bed
printf 'chr1\t16380\t16383\te1\nchr1\t16385\t16390\te1\nchr1\t131070\t131071\te2\nchr1\t131073\t131080\te2\nchr1\t1048570\t1048575\te3\nchr1\t1048577\t1048580\te3\nchr1\t16385\t16390\te4\nchr1\t131060\t131071\te5\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t019: B in a coarse bin and B over a 64 Mb edge
#  B chr2:0-600000000 crosses 2^29, so it is stored in a 2^32
#  bp bin, the second-coarsest level; it removes big1 and big2
#  and only touches big4. B
#  chr3:67108860-67108870 straddles 2^26 and cuts big3 in two.
###########################################################
echo -e "    subtract.new.t019...\c"
printf 'chr2\t300000000\t300000100\tbig1\nchr2\t1000\t1010\tbig2\nchr2\t600000000\t600000010\tbig4\nchr3\t67108800\t67108900\tbig3\n' > a.tmp.bed
printf 'chr2\t0\t600000000\nchr3\t67108860\t67108870\n' > b.tmp.bed
printf 'chr2\t600000000\t600000010\tbig4\nchr3\t67108800\t67108860\tbig3\nchr3\t67108870\t67108900\tbig3\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t020: coordinates above 2^31
#  h1 only touches B 3000000100-3000000200 and is printed
#  unchanged; h2 lies inside B 3999999000-4000000100 and is
#  removed. h3 lies inside B 4294967200-4294967400, which
#  crosses 2^32 and so is stored in the coarsest bin level;
#  h3 is removed.
###########################################################
echo -e "    subtract.new.t020...\c"
printf 'chr1\t3000000000\t3000000100\th1\nchr1\t4000000000\t4000000050\th2\nchr1\t4294967290\t4294967300\th3\n' > a.tmp.bed
printf 'chr1\t3000000100\t3000000200\nchr1\t3999999000\t4000000100\nchr1\t4294967200\t4294967400\n' > b.tmp.bed
printf 'chr1\t3000000000\t3000000100\th1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t021: zero-length A records
#  z1 (100-100) lies inside B 50-150 and is removed. z2
#  (900-900) has no B nearby and is printed unchanged. z3
#  (500-500) has an identical zero-length B and is removed.
###########################################################
echo -e "    subtract.new.t021...\c"
printf 'chr1\t100\t100\tz1\nchr1\t900\t900\tz2\nchr1\t500\t500\tz3\n' > a.tmp.bed
printf 'chr1\t50\t150\nchr1\t500\t500\n' > b.tmp.bed
printf 'chr1\t900\t900\tz2\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t023: BED5 A keeps name and score on each piece
###########################################################
echo -e "    subtract.new.t023...\c"
printf 'chr1\t100\t200\tn1\t9\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t150\tn1\t9\nchr1\t160\t200\tn1\t9\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t024: BED6 plus extra columns keeps every column
###########################################################
echo -e "    subtract.new.t024...\c"
printf 'chr1\t100\t200\tn1\t9\t+\textra1\textra2\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t150\tn1\t9\t+\textra1\textra2\nchr1\t160\t200\tn1\t9\t+\textra1\textra2\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t025: bedGraph A keeps its value on each piece
###########################################################
echo -e "    subtract.new.t025...\c"
printf 'chr1\t100\t200\t0.75\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t150\t0.75\nchr1\t160\t200\t0.75\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t026: BED12 A untouched or fully removed
#  tx1 has no overlapping B and is printed with all 12 fields;
#  tx2 lies inside B and is removed.
###########################################################
echo -e "    subtract.new.t026...\c"
printf 'chr1\t100\t200\ttx1\t0\t+\t110\t190\t255,0,0\t2\t20,20,\t0,80,\nchr1\t500\t600\ttx2\t0\t-\t500\t600\t0\t2\t10,10,\t0,90,\n' > a.tmp.bed
printf 'chr1\t400\t700\n' > b.tmp.bed
printf 'chr1\t100\t200\ttx1\t0\t+\t110\t190\t255,0,0\t2\t20,20,\t0,80,\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t027: GFF A pieces are printed 1-based
#  g1 is 101-200 (0-based 100-200); B 120-130 removes 1-based
#  121-130, leaving 101-120 and 131-200. g5 is 125-140
#  (0-based 124-140) and loses 1-based 125-130: 131-140.
###########################################################
echo -e "    subtract.new.t027...\c"
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\nchr1\tsrc\texon\t125\t140\t.\t-\t.\tID=g5\n' > a.tmp.gff
printf 'chr1\t120\t130\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t101\t120\t.\t+\t.\tID=g1\nchr1\tsrc\tgene\t131\t200\t.\t+\t.\tID=g1\nchr1\tsrc\texon\t131\t140\t.\t-\t.\tID=g5\n' > exp
$BT subtract -a a.tmp.gff -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  subtract.new.t028: GFF A with -A at 1-based edges
#  B 120-130 covers 1-based 121-130. g1 (101-200) overlaps it;
#  g2 (301-400) does not; g3 (131-140) only touches it; g4
#  (130-130) shares base 130. -A removes g1 and g4.
###########################################################
echo -e "    subtract.new.t028...\c"
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\nchr1\tsrc\tgene\t301\t400\t.\t+\t.\tID=g2\nchr1\tsrc\tgene\t131\t140\t.\t+\t.\tID=g3\nchr1\tsrc\tgene\t130\t130\t.\t+\t.\tID=g4\n' > a.tmp.gff
printf 'chr1\t120\t130\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t301\t400\t.\t+\t.\tID=g2\nchr1\tsrc\tgene\t131\t140\t.\t+\t.\tID=g3\n' > exp
$BT subtract -a a.tmp.gff -b b.tmp.bed -A > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  subtract.new.t029: GFF B is converted from 1-based
#  GFF 101-200 is 0-based 100-200. A 90-110 keeps 90-100;
#  A 200-210 only touches it and is kept whole.
###########################################################
echo -e "    subtract.new.t029...\c"
printf 'chr1\t90\t110\nchr1\t200\t210\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n' > b.tmp.gff
printf 'chr1\t90\t100\nchr1\t200\t210\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.gff > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  subtract.new.t030: VCF A records covered by B are removed
#  rs1 (POS 125, 0-based 124-125) is inside B 120-130; rs2
#  (POS 150, REF ACGT, 0-based 149-153) is inside B 140-160;
#  rs3 (POS 300) is not covered. The header is not printed.
###########################################################
echo -e "    subtract.new.t030...\c"
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t125\trs1\tA\tG\t50\tPASS\tDP=10\nchr1\t150\trs2\tACGT\tA\t50\tPASS\tDP=5\nchr1\t300\trs3\tC\tT\t50\tPASS\tDP=5\n' > a.tmp.vcf
printf 'chr1\t120\t130\nchr1\t140\t160\n' > b.tmp.bed
printf 'chr1\t300\trs3\tC\tT\t50\tPASS\tDP=5\n' > exp
$BT subtract -a a.tmp.vcf -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  subtract.new.t031: VCF A with -header keeps the VCF header
#  Same input as t030; both header lines precede rs3.
###########################################################
echo -e "    subtract.new.t031...\c"
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t125\trs1\tA\tG\t50\tPASS\tDP=10\nchr1\t150\trs2\tACGT\tA\t50\tPASS\tDP=5\nchr1\t300\trs3\tC\tT\t50\tPASS\tDP=5\n' > a.tmp.vcf
printf 'chr1\t120\t130\nchr1\t140\t160\n' > b.tmp.bed
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t300\trs3\tC\tT\t50\tPASS\tDP=5\n' > exp
$BT subtract -a a.tmp.vcf -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  subtract.new.t032: VCF B spans POS-1 to POS-1+len(REF)
#  rs1 removes base 124, rs2 (REF ACGT) removes 149-152, rs3
#  lies outside A. A 120-160 -> 120-124, 125-149, 153-160.
###########################################################
echo -e "    subtract.new.t032...\c"
printf 'chr1\t120\t160\n' > a.tmp.bed
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t125\trs1\tA\tG\t50\tPASS\tDP=10\nchr1\t150\trs2\tACGT\tA\t50\tPASS\tDP=5\nchr1\t300\trs3\tC\tT\t50\tPASS\tDP=5\n' > b.tmp.vcf
printf 'chr1\t120\t124\nchr1\t125\t149\nchr1\t153\t160\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.vcf > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  subtract.new.t033: gzipped A and B read like plain text
#  Same records as t024 (BED6+2 A cut by 150-160).
###########################################################
echo -e "    subtract.new.t033...\c"
printf 'chr1\t100\t200\tn1\t9\t+\textra1\textra2\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
gzip -c a.tmp.bed > a.tmp.bed.gz
gzip -c b.tmp.bed > b.tmp.bed.gz
printf 'chr1\t100\t150\tn1\t9\t+\textra1\textra2\nchr1\t160\t200\tn1\t9\t+\textra1\textra2\n' > exp
$BT subtract -a a.tmp.bed.gz -b b.tmp.bed.gz > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed a.tmp.bed.gz b.tmp.bed.gz

###########################################################
#  subtract.new.t034: BAM A gives BAM output
#  r1 (100-150) is partly covered by B 120-130, r2 (200-230)
#  lies inside B 190-240, r3 (400-420) has no overlap. r2 is
#  removed; r1 and r3 are written. Pins the current output: a
#  partly covered read is written whole and once, since BAM
#  records are not trimmed.
###########################################################
echo -e "    subtract.new.t034...\c"
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*\nr2\t16\tchr1\t201\t30\t30M\t*\t0\t0\t*\t*\nr3\t0\tchr1\t401\t60\t20M\t*\t0\t0\t*\t*\n' > a.tmp.sam
$htsutil samtobam a.tmp.sam a.tmp.bam
printf 'chr1\t120\t130\nchr1\t190\t240\n' > b.tmp.bed
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*\nr3\t0\tchr1\t401\t60\t20M\t*\t0\t0\t*\t*\n' > exp
$BT subtract -a a.tmp.bam -b b.tmp.bed > o.tmp.bam
$htsutil viewbam o.tmp.bam > obs
check obs exp
rm obs exp a.tmp.sam a.tmp.bam o.tmp.bam b.tmp.bed

###########################################################
#  subtract.new.t035: BAM A with -A drops every overlapped read
#  Same input as t034; only r3 has no overlap.
###########################################################
echo -e "    subtract.new.t035...\c"
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*\nr2\t16\tchr1\t201\t30\t30M\t*\t0\t0\t*\t*\nr3\t0\tchr1\t401\t60\t20M\t*\t0\t0\t*\t*\n' > a.tmp.sam
$htsutil samtobam a.tmp.sam a.tmp.bam
printf 'chr1\t120\t130\nchr1\t190\t240\n' > b.tmp.bed
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr3\t0\tchr1\t401\t60\t20M\t*\t0\t0\t*\t*\n' > exp
$BT subtract -a a.tmp.bam -b b.tmp.bed -A > o.tmp.bam
$htsutil viewbam o.tmp.bam > obs
check obs exp
rm obs exp a.tmp.sam a.tmp.bam o.tmp.bam b.tmp.bed

###########################################################
#  subtract.new.t036: BAM B without -split uses the read span
#  r1 is 50M at 100-150; r2 is 20M30N20M at 200-270. A 90-110
#  keeps 90-100; A 225-245 lies in r2's span and is removed.
###########################################################
echo -e "    subtract.new.t036...\c"
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*\nr2\t16\tchr1\t201\t30\t20M30N20M\t*\t0\t0\t*\t*\n' > b.tmp.sam
$htsutil samtobam b.tmp.sam b.tmp.bam
printf 'chr1\t90\t110\nchr1\t225\t245\n' > a.tmp.bed
printf 'chr1\t90\t100\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bam > obs
check obs exp
rm obs exp b.tmp.sam b.tmp.bam a.tmp.bed

###########################################################
#  subtract.new.t037: BAM B with -split uses aligned blocks
#  Same input as t036. A 225-245 lies in r2's 30N gap
#  (220-250), so with -split it has no overlap and is kept.
###########################################################
echo -e "    subtract.new.t037...\c"
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*\nr2\t16\tchr1\t201\t30\t20M30N20M\t*\t0\t0\t*\t*\n' > b.tmp.sam
$htsutil samtobam b.tmp.sam b.tmp.bam
printf 'chr1\t90\t110\nchr1\t225\t245\n' > a.tmp.bed
printf 'chr1\t90\t100\nchr1\t225\t245\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bam -split > obs
check obs exp
rm obs exp b.tmp.sam b.tmp.bam a.tmp.bed

###########################################################
#  subtract.new.t038: -f is tested against each B separately
#  B 100-130 and 150-180 each cover 30% of A (60% together).
#  Neither reaches -f 0.5, so A is kept whole.
###########################################################
echo -e "    subtract.new.t038...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t130\nchr1\t150\t180\n' > b.tmp.bed
printf 'chr1\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t039: -f includes a B at exactly the fraction
#  B 100-150 covers exactly 0.5 of A and is removed (the test
#  is >=); B 180-190 covers 0.1 and is ignored: 150-200.
###########################################################
echo -e "    subtract.new.t039...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t150\nchr1\t180\t190\n' > b.tmp.bed
printf 'chr1\t150\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t040: -f 1.0 needs one B covering all of A
#  A1 has an identical B and is removed; A2's B covers 99 of
#  its 100 bp and is ignored; A3's B contains it.
###########################################################
echo -e "    subtract.new.t040...\c"
printf 'chr1\t100\t200\nchr1\t300\t400\nchr1\t500\t600\n' > a.tmp.bed
printf 'chr1\t100\t200\nchr1\t300\t399\nchr1\t450\t650\n' > b.tmp.bed
printf 'chr1\t300\t400\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 1.0 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t041: -F is a fraction of B
#  B1 120-130 lies wholly in A (fraction 1.0); B2 150-400 has
#  50 of 250 bp in A (0.2). With -F 0.5 only B1 is removed.
###########################################################
echo -e "    subtract.new.t041...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t120\t130\nchr1\t150\t400\n' > b.tmp.bed
printf 'chr1\t100\t120\nchr1\t130\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -F 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t042: -F 1.0 keeps only B wholly inside A
#  B1 120-130 and B3 160-200 are inside A; B2 190-210 has 10
#  of 20 bp in A. Removing B1 and B3 leaves 100-120, 130-160.
###########################################################
echo -e "    subtract.new.t042...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t120\t130\nchr1\t190\t210\nchr1\t160\t200\n' > b.tmp.bed
printf 'chr1\t100\t120\nchr1\t130\t160\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -F 1.0 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t043: -f with and without -r
#  A 100-200. B1 150-400: 0.5 of A, 0.2 of B. B2 90-110: 0.1 of
#  A. B3 175-215: 0.25 of A, 0.625 of B.
#  -f 0.25 -r keeps only B3 -> 100-175.
#  -f 0.25 alone keeps B1 and B3 -> 100-150.
###########################################################
echo -e "    subtract.new.t043...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t400\nchr1\t90\t110\nchr1\t175\t215\n' > b.tmp.bed
printf 'chr1\t100\t175\nchr1\t100\t150\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 0.25 -r > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 0.25 >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t044: -f with -F needs both; -e needs either
#  A 100-200. B1 150-400 covers 0.5 of A but only 0.2 of B1.
#  B2 90-110 covers 0.1 of A and 0.5 of B2.
#  -f 0.5 -F 0.5: neither qualifies -> A whole.
#  -f 0.5 -F 0.5 -e: both qualify -> 110-150.
###########################################################
echo -e "    subtract.new.t044...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t400\nchr1\t90\t110\n' > b.tmp.bed
printf 'chr1\t100\t200\nchr1\t110\t150\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.5 > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.5 -e >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t045: -A removes overlapped A and keeps the rest
#  k1 has no B; k2 shares base 390 with B; k3 only touches
#  B chr2:200-300. Output keeps k1 and k3 in input order.
###########################################################
echo -e "    subtract.new.t045...\c"
printf 'chr1\t100\t200\tk1\nchr1\t300\t400\tk2\nchr2\t100\t200\tk3\n' > a.tmp.bed
printf 'chr1\t390\t391\nchr2\t200\t300\nchr1\t500\t600\n' > b.tmp.bed
printf 'chr1\t100\t200\tk1\nchr2\t100\t200\tk3\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -A > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t046: -A with a zero-length B inside A
#  The zero-length B at 150 overlaps A 100-200 and removes it;
#  A 300-400 is kept.
###########################################################
echo -e "    subtract.new.t046...\c"
printf 'chr1\t100\t200\nchr1\t300\t400\n' > a.tmp.bed
printf 'chr1\t150\t150\n' > b.tmp.bed
printf 'chr1\t300\t400\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -A > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t047: -A with -F
#  A1 100-200 overlaps B1 150-400 with 0.2 of B1; A2
#  1000-1100 overlaps B2 1050-1110 with 50/60 of B2. With
#  -F 0.5 only A2 is removed.
###########################################################
echo -e "    subtract.new.t047...\c"
printf 'chr1\t100\t200\nchr1\t1000\t1100\n' > a.tmp.bed
printf 'chr1\t150\t400\nchr1\t1050\t1110\n' > b.tmp.bed
printf 'chr1\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -A -F 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t048: -A with -f -r
#  A1's B 150-400 covers 0.5 of A1 but A1 is only 0.2 of it,
#  so -f 0.5 -r keeps A1. A2's B 1040-1110 covers 0.6 of A2
#  and A2 is 60/70 of it, so A2 is removed.
###########################################################
echo -e "    subtract.new.t048...\c"
printf 'chr1\t100\t200\nchr1\t1000\t1100\n' > a.tmp.bed
printf 'chr1\t150\t400\nchr1\t1040\t1110\n' > b.tmp.bed
printf 'chr1\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -A -f 0.5 -r > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t049: -A with -s and with -S
#  A records on +, - and unknown (.) strands; B is on +.
#  -s removes only p; -S removes only m; d is never removed
#  because an unknown strand matches neither rule.
###########################################################
echo -e "    subtract.new.t049...\c"
printf 'chr1\t100\t200\tp\t0\t+\nchr1\t100\t200\tm\t0\t-\nchr1\t100\t200\td\t0\t.\n' > a.tmp.bed
printf 'chr1\t150\t160\tb\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\tm\t0\t-\nchr1\t100\t200\td\t0\t.\nchr1\t100\t200\tp\t0\t+\nchr1\t100\t200\td\t0\t.\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -A -s > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -A -S >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t050: -A with several -b files
#  b1 overlaps the first A record and b2 the third; only
#  chr1:300-400 is kept.
###########################################################
echo -e "    subtract.new.t050...\c"
printf 'chr1\t100\t200\nchr1\t300\t400\nchr1\t500\t600\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b1.tmp.bed
printf 'chr1\t550\t560\n' > b2.tmp.bed
printf 'chr1\t300\t400\n' > exp
$BT subtract -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed -A > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  subtract.new.t051: BED3 A in the intron of a BED12 B
#  B spans 100-200 with blocks 100-120 and 180-200. Without
#  -split both A records lie inside B's span and are removed.
#  With -split i1 (130-150, in the intron) is kept and i2
#  (110-115, in a block) is removed.
###########################################################
echo -e "    subtract.new.t051...\c"
printf 'chr1\t130\t150\ti1\nchr1\t110\t115\ti2\n' > a.tmp.bed
printf 'chr1\t100\t200\ttx\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\n' > b.tmp.bed
printf 'chr1\t130\t150\ti1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -split >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t052: BED12 A with B only in its intron
#  A has blocks 100-120 and 180-200; B 130-150 touches neither
#  block, so with -split A is printed unchanged, and -split -A
#  keeps it too.
###########################################################
echo -e "    subtract.new.t052...\c"
printf 'chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\n' > a.tmp.bed
printf 'chr1\t130\t150\n' > b.tmp.bed
printf 'chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\nchr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -split > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -split -A >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t053: -split measures -f against A's blocks
#  A blocks 100-120 and 180-200 total 40 bp; B 110-130 covers
#  10 of them (0.25) but 20 of A's 100 bp span (0.2).
#  -A -f 0.25 keeps A; -A -f 0.25 -split removes it.
###########################################################
echo -e "    subtract.new.t053...\c"
printf 'chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\n' > a.tmp.bed
printf 'chr1\t110\t130\n' > b.tmp.bed
printf 'chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t20,20,\t0,80,\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -A -f 0.25 > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -A -f 0.25 -split >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t054: -N totals all B where -A tests each B
#  B 100-130 and 120-160 cover 30% and 40% of A. Their total
#  is above 0.5 whether summed (70%) or merged (60%), so
#  -N -f 0.5 removes A; -A -f 0.5 keeps it.
###########################################################
echo -e "    subtract.new.t054...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t130\nchr1\t120\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -N -f 0.5 > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -A -f 0.5 >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t056: -N with -s sums only same-strand B
#  B1 100-140 is on A's strand, B2 140-180 on the other.
#  -N -f 0.5 -s sees 40% and keeps A; -N -f 0.5 sees 80%
#  and removes it.
###########################################################
echo -e "    subtract.new.t056...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf 'chr1\t100\t140\tb\t0\t+\nchr1\t140\t180\tc\t0\t-\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\t0\t+\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -N -f 0.5 -s > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -N -f 0.5 >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t057: -N with -F filters B before summing
#  B1 100-130 lies in A; B2 140-400 has 60 of 260 bp in A.
#  -N -f 0.5 sums 90% and removes A. Adding -F 0.5 drops B2,
#  leaving 30%, so A is kept.
###########################################################
echo -e "    subtract.new.t057...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t130\nchr1\t140\t400\n' > b.tmp.bed
printf 'chr1\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -N -f 0.5 > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -N -f 0.5 -F 0.5 >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t058: -N keeps A records with no overlap
#  chr1:100-200 has no B and is kept; chr2:100-200 is fully
#  covered (1.0 > 0.5) and removed.
###########################################################
echo -e "    subtract.new.t058...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr2\t100\t200\n' > b.tmp.bed
printf 'chr1\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -N -f 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t059: -s subtracts only same-strand B
#  A on +, - and unknown (.). B: bp 110-120 (+), bm 130-140
#  (-), bd 150-160 (.). p loses 110-120, m loses 130-140, d
#  loses nothing: an unknown strand never matches, as stated
#  in Record::sameChromIntersects.
###########################################################
echo -e "    subtract.new.t059...\c"
printf 'chr1\t100\t200\tp\t0\t+\nchr1\t100\t200\tm\t0\t-\nchr1\t100\t200\td\t0\t.\n' > a.tmp.bed
printf 'chr1\t110\t120\tbp\t0\t+\nchr1\t130\t140\tbm\t0\t-\nchr1\t150\t160\tbd\t0\t.\n' > b.tmp.bed
printf 'chr1\t100\t110\tp\t0\t+\nchr1\t120\t200\tp\t0\t+\nchr1\t100\t130\tm\t0\t-\nchr1\t140\t200\tm\t0\t-\nchr1\t100\t200\td\t0\t.\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -s > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t060: -S subtracts only opposite-strand B
#  Same input as t059: p loses 130-140, m loses 110-120, d
#  loses nothing.
###########################################################
echo -e "    subtract.new.t060...\c"
printf 'chr1\t100\t200\tp\t0\t+\nchr1\t100\t200\tm\t0\t-\nchr1\t100\t200\td\t0\t.\n' > a.tmp.bed
printf 'chr1\t110\t120\tbp\t0\t+\nchr1\t130\t140\tbm\t0\t-\nchr1\t150\t160\tbd\t0\t.\n' > b.tmp.bed
printf 'chr1\t100\t130\tp\t0\t+\nchr1\t140\t200\tp\t0\t+\nchr1\t100\t110\tm\t0\t-\nchr1\t120\t200\tm\t0\t-\nchr1\t100\t200\td\t0\t.\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -S > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t061: -s with GFF strands
#  Both GFF records span 101-200; B (0-based 120-130) is on
#  -. The + record is kept whole; the - record loses 1-based
#  121-130.
###########################################################
echo -e "    subtract.new.t061...\c"
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=p\nchr1\tsrc\tgene\t101\t200\t.\t-\t.\tID=m\n' > a.tmp.gff
printf 'chr1\t120\t130\tb\t0\t-\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=p\nchr1\tsrc\tgene\t101\t120\t.\t-\t.\tID=m\nchr1\tsrc\tgene\t131\t200\t.\t-\t.\tID=m\n' > exp
$BT subtract -a a.tmp.gff -b b.tmp.bed -s > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  subtract.new.t062: -s with a strandless (BED5) B is an error
###########################################################
echo -e "    subtract.new.t062...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf 'chr1\t150\t160\tb\t5\n' > b.tmp.bed
printf '***** ERROR: stranded subtract requested, but input file b.tmp.bed does not have strands. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -s > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t063: -s with VCF input is an error
###########################################################
echo -e "    subtract.new.t063...\c"
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t125\trs1\tA\tG\t50\tPASS\tDP=10\n' > a.tmp.vcf
printf 'chr1\t120\t130\tb\t0\t+\n' > b.tmp.bed
printf '***** ERROR: stranded subtract not supported for VCF file a.tmp.vcf. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.vcf -b b.tmp.bed -s > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.vcf b.tmp.bed

###########################################################
#  subtract.new.t064: -s and -S together are an error
###########################################################
echo -e "    subtract.new.t064...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf 'chr1\t150\t160\tb\t0\t+\n' > b.tmp.bed
printf '***** ERROR: request -s for sameStrand, or -S for diffStrand, not both. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -s -S > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t065: three -b files whose B overlap each other
#  b1 110-120, 115-125; b2 150-160, chr2 190-250; b3 118-152.
#  On chr1 the union is 110-160, leaving 100-110 and 160-200;
#  chr2 keeps 100-190.
###########################################################
echo -e "    subtract.new.t065...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t110\t120\nchr1\t115\t125\n' > b1.tmp.bed
printf 'chr1\t150\t160\nchr2\t190\t250\n' > b2.tmp.bed
printf 'chr1\t118\t152\n' > b3.tmp.bed
printf 'chr1\t100\t110\nchr1\t160\t200\nchr2\t100\t190\n' > exp
$BT subtract -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed b3.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed b3.tmp.bed

###########################################################
#  subtract.new.t066: repeated -b flags act like one list
#  -b b1 -b b2 with the b1/b2 of t065: chr1 loses 110-125 and
#  150-160; chr2 loses 190-200.
###########################################################
echo -e "    subtract.new.t066...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t110\t120\nchr1\t115\t125\n' > b1.tmp.bed
printf 'chr1\t150\t160\nchr2\t190\t250\n' > b2.tmp.bed
printf 'chr1\t100\t110\nchr1\t125\t150\nchr1\t160\t200\nchr2\t100\t190\n' > exp
$BT subtract -a a.tmp.bed -b b1.tmp.bed -b b2.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  subtract.new.t067: -f with several -b files tests each B
#  b1 100-130 covers 0.3 of A and b2 130-200 covers 0.7.
#  With -f 0.5 only b2 is removed: 100-130 remains.
###########################################################
echo -e "    subtract.new.t067...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t130\n' > b1.tmp.bed
printf 'chr1\t130\t200\n' > b2.tmp.bed
printf 'chr1\t100\t130\n' > exp
$BT subtract -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed -f 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  subtract.new.t068: -b files of different formats
#  b1 is BED (100-110); b2 is GFF 191-200 (0-based 190-200).
#  A 100-200 keeps 110-190.
###########################################################
echo -e "    subtract.new.t068...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t100\t110\n' > b1.tmp.bed
printf 'chr1\tsrc\tgene\t191\t200\t.\t+\t.\tID=g\n' > b2.tmp.gff
printf 'chr1\t110\t190\n' > exp
$BT subtract -a a.tmp.bed -b b1.tmp.bed b2.tmp.gff > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.gff

###########################################################
#  subtract.new.t069: -sorted on several chromosomes
#  chr1 and chr3 have no B; chr4 is only in B; chr2 loses
#  150-160.
###########################################################
echo -e "    subtract.new.t069...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\nchr3\t100\t200\n' > a.tmp.bed
printf 'chr2\t150\t160\nchr4\t1\t5\n' > b.tmp.bed
printf 'chr1\t100\t200\nchr2\t100\t150\nchr2\t160\t200\nchr3\t100\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -sorted > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t070: -sorted -g with chr2 ordered before chr1
###########################################################
echo -e "    subtract.new.t070...\c"
printf 'chr2\t100\t200\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr2\t150\t160\nchr1\t120\t130\n' > b.tmp.bed
printf 'chr2\t1000\nchr1\t1000\n' > g.tmp.txt
printf 'chr2\t100\t150\nchr2\t160\t200\nchr1\t100\t120\nchr1\t130\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -sorted -g g.tmp.txt > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed g.tmp.txt

###########################################################
#  subtract.new.t071: -sorted with an unsorted A is an error
#  The message names the file and the first out-of-order record.
###########################################################
echo -e "    subtract.new.t071...\c"
printf 'chr1\t300\t400\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'Error: Sorted input specified, but the file a.tmp.bed has the following out of order record\nchr1\t100\t200\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -sorted > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t072: -sorted -g rejects order unlike the genome
#  The genome file puts chrX first but B lists chr1 first.
#  Pins the current output: B's record is the one reported.
###########################################################
echo -e "    subtract.new.t072...\c"
printf 'chr1\t100\t200\nchrX\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\nchrX\t150\t160\n' > b.tmp.bed
printf 'chrX\t1000\nchr1\t1000\n' > g.tmp.txt
printf 'Error: Sorted input specified, but the file b.tmp.bed has the following record with a different sort order than the genomeFile g.tmp.txt\nchrX\t150\t160\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -sorted -g g.tmp.txt > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed g.tmp.txt

###########################################################
#  subtract.new.t073: -sorted with two -b files
#  b1 110-120 and b2 150-160 cut chr1 into three pieces; b2
#  chr2:190-210 trims chr2 to 100-190.
###########################################################
echo -e "    subtract.new.t073...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t110\t120\n' > b1.tmp.bed
printf 'chr1\t150\t160\nchr2\t190\t210\n' > b2.tmp.bed
printf 'chr1\t100\t110\nchr1\t120\t150\nchr1\t160\t200\nchr2\t100\t190\n' > exp
$BT subtract -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed -sorted > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  subtract.new.t074: -g without -sorted leaves output unchanged
#  Unsorted A (chr2 before chr1) with a genome file listing
#  both chromosomes gives the ordinary unsorted result.
###########################################################
echo -e "    subtract.new.t074...\c"
printf 'chr2\t100\t200\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t1000\nchr2\t1000\n' > g.tmp.txt
printf 'chr2\t100\t200\nchr1\t100\t150\nchr1\t160\t200\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -g g.tmp.txt > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed g.tmp.txt

###########################################################
#  subtract.new.t075: chr1 vs 1 naming warns on stderr
#  The names differ, so nothing overlaps and A is kept. Pins
#  the current output: B is loaded before A is read, so the
#  warning names A's record, and it is printed when first
#  tripped and again at exit.
###########################################################
echo -e "    subtract.new.t075...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf '1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\n--- stderr\n***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t100\t200\n\n***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t100\t200\n\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs 2> err.tmp
echo "--- stderr" >> obs
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t076: -nonamecheck silences that warning
#  Same input as t075, unsorted: stdout as before, stderr empty.
###########################################################
echo -e "    subtract.new.t076...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf '1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\n--- stderr\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -nonamecheck > obs 2> err.tmp
echo "--- stderr" >> obs
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t077: -sorted -nonamecheck with chr01 vs chr1
#  No overlap (different names), no warning.
###########################################################
echo -e "    subtract.new.t077...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr01\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\n--- stderr\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -sorted -nonamecheck > obs 2> err.tmp
echo "--- stderr" >> obs
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t078: -header prints A's header, not B's
#  Without -header the '#' and 'track' lines are dropped; with
#  it they come first, once. B's own header is never printed.
###########################################################
echo -e "    subtract.new.t078...\c"
printf '#comment line\ntrack name=foo\nchr1\t100\t200\tr1\nchr2\t10\t20\tr2\n' > a.tmp.bed
printf '#b header\nchr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t150\tr1\nchr1\t160\t200\tr1\nchr2\t10\t20\tr2\n#comment line\ntrack name=foo\nchr1\t100\t150\tr1\nchr1\t160\t200\tr1\nchr2\t10\t20\tr2\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
$BT subtract -a a.tmp.bed -b b.tmp.bed -header >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t079: -header is printed when every A is removed
###########################################################
echo -e "    subtract.new.t079...\c"
printf '#h\nchr1\t150\t160\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '#h\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t080: -v is rejected by subtract
###########################################################
echo -e "    subtract.new.t080...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: -v option is not valid for subtract. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -v > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t081: -f outside (0, 1] is an error
#  Checked for -f 0 and -f 1.5.
###########################################################
echo -e "    subtract.new.t081...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: -f must be in the range (0.0, 1.0]. *****\nexit 1\n***** ERROR: -f must be in the range (0.0, 1.0]. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 0 > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 1.5 >> out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t082: -F outside (0, 1] is an error
###########################################################
echo -e "    subtract.new.t082...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: -F must be in the range (0.0, 1.0]. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -F 0 > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t083: -r with -F is an error
###########################################################
echo -e "    subtract.new.t083...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: -r must be used solely with -f. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.5 -r > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t084: -e without -f is an error
###########################################################
echo -e "    subtract.new.t084...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: -e must be used with -f. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -F 0.5 -e > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t085: -a without -b is an error
###########################################################
echo -e "    subtract.new.t085...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf '***** ERROR: query and database files not specified. *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed

###########################################################
#  subtract.new.t086: a missing input file is an error
#  Checked for a missing A and a missing B.
###########################################################
echo -e "    subtract.new.t086...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'Error: Unable to open file missing.tmp.bed. Exiting.\nexit 1\nError: Unable to open file missing.tmp.bed. Exiting.\nexit 1\n' > exp
$BT subtract -a missing.tmp.bed -b a.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
$BT subtract -a a.tmp.bed -b missing.tmp.bed >> obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  subtract.new.t087: an unknown option is an error
###########################################################
echo -e "    subtract.new.t087...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf '***** ERROR: Unrecognized parameter: -bogus *****\nexit 1\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed -bogus > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
grep -e 'ERROR' -e '^exit' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  subtract.new.t088: empty B keeps A; empty A prints nothing
#  Both runs exit 0.
###########################################################
echo -e "    subtract.new.t088...\c"
printf 'chr1\t100\t200\nchr2\t10\t20\n' > a.tmp.bed
printf '' > e.tmp.bed
printf 'chr1\t100\t200\nchr2\t10\t20\nexit 0\nexit 0\n' > exp
$BT subtract -a a.tmp.bed -b e.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
$BT subtract -a e.tmp.bed -b a.tmp.bed >> obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed e.tmp.bed

###########################################################
#  subtract.new.t089: BAM A with -bed prints reads as BED12
#  r2 (200-230) lies inside B 190-240 and is removed; r3
#  (20M at 400-420) has no overlap and is printed in the
#  BED12 form of a BAM record.
###########################################################
echo -e "    subtract.new.t089...\c"
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr2\t16\tchr1\t201\t30\t30M\t*\t0\t0\t*\t*\nr3\t0\tchr1\t401\t60\t20M\t*\t0\t0\t*\t*\n' > a.tmp.sam
$htsutil samtobam a.tmp.sam a.tmp.bam
printf 'chr1\t190\t240\n' > b.tmp.bed
printf 'chr1\t400\t420\tr3\t60\t+\t400\t420\t0,0,0\t1\t20,\t0,\n' > exp
$BT subtract -a a.tmp.bam -b b.tmp.bed -bed > obs
check obs exp
rm obs exp a.tmp.sam a.tmp.bam b.tmp.bed

###########################################################
#  subtract.new.t090: a VCF structural variant as B
#  sv2 <DUP> at POS 201 with END=300 covers far more than its
#  1 bp REF. A 210-260 lies inside it under any reading of
#  its first and last base and is removed; A 100-150 is far
#  from it and kept. (The SVLEN/END span disagreement is in
#  bugs.md and is not asserted.)
###########################################################
echo -e "    subtract.new.t090...\c"
printf 'chr1\t100\t150\nchr1\t210\t260\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t201\tsv2\tN\t<DUP>\t50\tPASS\tSVTYPE=DUP;END=300\n' > b.tmp.vcf
printf 'chr1\t100\t150\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.vcf > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  subtract.new.t091: A on either side of the first fine-bin
#  edge (16384), each needing B from both fine bins
#  f1 16383-16390 starts on the last base of fine bin 0; f2
#  16380-16385 ends on the first base of fine bin 1. B
#  16383-16384 sits in bin 0 and B 16384-16385 in bin 1, and
#  each A overlaps both. f1 keeps 16385-16390, f2 keeps
#  16380-16383. The B are kept apart from t018's, whose
#  16383-16385 would cover both bases on its own.
###########################################################
echo -e "    subtract.new.t091...\c"
printf 'chr1\t16383\t16390\tf1\nchr1\t16380\t16385\tf2\n' > a.tmp.bed
printf 'chr1\t16383\t16384\nchr1\t16384\t16385\n' > b.tmp.bed
printf 'chr1\t16385\t16390\tf1\nchr1\t16380\t16383\tf2\n' > exp
$BT subtract -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

[[ $FAILURES -eq 0 ]] || exit 1;
