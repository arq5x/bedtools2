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

###########################################################
#  Test a basic self intersection
###########################################################
echo -e "    jaccard.t01...\c"
echo \
"intersection	union	jaccard	n_intersections
110	110	1	2" > exp
$BT jaccard -a a.bed -b a.bed > obs
check obs exp
rm obs exp


echo -e "    jaccard.t02...\c"
echo \
"intersection	union	jaccard	n_intersections
10	140	0.0714286	1" > exp
$BT jaccard -a a.bed -b b.bed > obs
check obs exp
rm obs exp

echo -e "    jaccard.t03...\c"
echo \
"intersection	union	jaccard	n_intersections
10	200	0.05	1" > exp
$BT jaccard -a a.bed -b c.bed > obs
check obs exp
rm obs exp

# TEST #4 IS DEPRECATED

###########################################################
#  Test stdin
###########################################################
echo -e "    jaccard.t05...\c"
echo \
"intersection	union	jaccard	n_intersections
10	140	0.0714286	1" > exp
cat a.bed | $BT jaccard -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test symmetry
###########################################################
echo -e "    jaccard.t06...\c"
$BT jaccard -a a.bed -b b.bed > obs1
$BT jaccard -a b.bed -b a.bed > obs2
check obs1 obs2
rm obs1 obs2

###########################################################
#  Test partially matching blocks without -split option.
###########################################################
echo -e "    jaccard.t07...\c"
echo \
"intersection	union	jaccard	n_intersections
10	50	0.2	1" > exp
$BT jaccard -a three_blocks_match.bed -b e.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test partially matching blocks with -split option.
###########################################################
echo -e "    jaccard.t08...\c"
echo \
"intersection	union	jaccard	n_intersections
5	35	0.142857	1" > exp
$BT jaccard -a three_blocks_match.bed -b e.bed -split > obs
check obs exp
rm obs exp

###########################################################
#  Test jaccard of Bam with Bam
###########################################################
echo -e "    jaccard.t09...\c"
echo \
"intersection	union	jaccard	n_intersections
10	150	0.0666667	1" > exp
$BT jaccard -a a.bam -b three_blocks_match.bam -bed > obs
check exp obs
rm exp obs

###########################################################
#  Test jaccard with mixed strand files
###########################################################
echo -e "    jaccard.t10...\c"
echo \
"intersection	union	jaccard	n_intersections
145	180	0.805556	2" >exp
$BT jaccard -a aMixedStrands.bed -b bMixedStrands.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test jaccard with mixed strand files, -s option
#  (match strand, either forward or reverse)
###########################################################
echo -e "    jaccard.t11...\c"
echo \
"intersection	union	jaccard	n_intersections
120	290	0.413793	4" >exp
$BT jaccard -a aMixedStrands.bed -b bMixedStrands.bed -s > obs
check obs exp
rm obs exp

###########################################################
#  Test jaccard with mixed strand files, -S + option
#  (match strand, forward only)
###########################################################
echo -e "    jaccard.t12...\c"
echo \
"intersection	union	jaccard	n_intersections
40	135	0.296296	2" >exp
$BT jaccard -a aMixedStrands.bed -b bMixedStrands.bed -S + > obs
check obs exp
rm obs exp

###########################################################
#  Test jaccard with mixed strand files, -S - option
#  (match strand, reverse only)
###########################################################
echo -e "    jaccard.t13...\c"
echo \
"intersection	union	jaccard	n_intersections
80	155	0.516129	2" > exp
$BT jaccard -a aMixedStrands.bed -b bMixedStrands.bed -S - > obs
check obs exp
rm obs exp

echo -e "    jaccard.t14...\c"
echo  "intersection	union	jaccard	n_intersections
1	3	0.333333	1" > exp

$BT jaccard -b a645.bed -a b645.bed > obs
check obs exp
rm obs 

echo -e "    jaccard.t15...\c"
$BT jaccard -a a645.bed -b b645.bed > obs
check obs exp
rm obs  exp

echo -e "    jaccard.t16...\c"
echo -e "intersection	union	jaccard	n_intersections
247800000	615800000	0.402403	4" > exp
$BT jaccard -a long.bed -b short.bed > obs
check obs exp
rm obs  exp

###########################################################
#
# NEW TESTS: output columns, merging, options, record types
#
# jaccard merges overlapping and book-ended intervals within
# each file, then reports: intersection = bases shared by the
# merged A and merged B intervals; union = total merged A
# length + total merged B length - intersection; jaccard =
# intersection / union, printed with 6 significant digits;
# n_intersections = number of merged A / merged B pairs that
# overlap.
#
###########################################################

###########################################################
#  jaccard.new.t001: disjoint intervals
#  A chr1:0-100, B chr1:200-300: intersection 0, union
#  100 + 100 - 0 = 200, jaccard 0, no intersections.
###########################################################
echo -e "    jaccard.new.t001...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t200\t300\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n0\t200\t0\t0\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t002: book-ended A and B do not intersect
#  A chr1:0-100, B chr1:100-200 share no base: 0 200 0 0.
###########################################################
echo -e "    jaccard.new.t002...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n0\t200\t0\t0\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t003: a partial overlap
#  A chr1:0-100, B chr1:50-150: intersection 50, union
#  100 + 100 - 50 = 150, jaccard 50/150 = 0.333333, 1 pair.
###########################################################
echo -e "    jaccard.new.t003...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n50\t150\t0.333333\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t004: the jaccard value is rounded to 6 significant
#  digits
#  A chr1:0-3, B chr1:1-3: intersection 2, union 3 + 2 - 2 = 3,
#  2/3 = 0.6666666... printed as 0.666667. pins the current output:
#  6 significant digits, rounded.
###########################################################
echo -e "    jaccard.new.t004...\c"
printf 'chr1\t0\t3\n' > a.tmp.bed
printf 'chr1\t1\t3\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n2\t3\t0.666667\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t005: overlapping A intervals are merged first
#  A chr1:0-100 and chr1:50-150 merge to chr1:0-150; B
#  chr1:0-150: intersection 150, union 150, jaccard 1, 1 pair.
###########################################################
echo -e "    jaccard.new.t005...\c"
printf 'chr1\t0\t100\nchr1\t50\t150\n' > a.tmp.bed
printf 'chr1\t0\t150\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n150\t150\t1\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t006: book-ended A intervals are merged too
#  A chr1:0-100 and chr1:100-200 merge to chr1:0-200; B
#  chr1:50-150: intersection 100, union 200 + 100 - 100 = 200,
#  jaccard 0.5, 1 pair. pins the current output: book-ended records
#  are merged before counting, so n_intersections is 1; the usage
#  text and jaccard.rst do not mention merging.
###########################################################
echo -e "    jaccard.new.t006...\c"
printf 'chr1\t0\t100\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t150\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n100\t200\t0.5\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t007: overlapping B intervals are merged
#  B chr1:0-50 and chr1:25-75 merge to chr1:0-75; A chr1:0-100:
#  intersection 75, union 100 + 75 - 75 = 100, jaccard 0.75.
###########################################################
echo -e "    jaccard.new.t007...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t0\t50\nchr1\t25\t75\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n75\t100\t0.75\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t008: one B covering several separate A intervals
#  A chr1:0-10, 20-30, 40-50 (30 bp) inside B chr1:0-100:
#  intersection 30, union 30 + 100 - 30 = 100, jaccard 0.3,
#  3 pairs.
###########################################################
echo -e "    jaccard.new.t008...\c"
printf 'chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\n' > a.tmp.bed
printf 'chr1\t0\t100\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n30\t100\t0.3\t3\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t009: one A overlapping two B intervals
#  A chr1:0-100; B chr1:10-20 and chr1:30-40: intersection 20,
#  union 100 + 20 - 20 = 100, jaccard 0.2, 2 pairs.
###########################################################
echo -e "    jaccard.new.t009...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\nchr1\t30\t40\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n20\t100\t0.2\t2\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t010: several chromosomes, including one only B has
#  A chr1:0-100 and chr2:0-100; B chr1:50-100, chr2:0-50 and
#  chr3:0-100. Intersection 50 + 50 = 100; union 200 + 200 - 100
#  = 300 (the chr3 B counts); jaccard 0.333333; 2 pairs.
###########################################################
echo -e "    jaccard.new.t010...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t100\nchr2\t0\t50\nchr3\t0\t100\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n100\t300\t0.333333\t2\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t011: A and B on different chromosomes
#  A chr1:0-100, B chr2:0-100: nothing shared, union 200.
###########################################################
echo -e "    jaccard.new.t011...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr2\t0\t100\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n0\t200\t0\t0\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t012: -f rejects an overlap below the fraction of A
#  A chr1:0-100, B chr1:50-150: 50 bp = 0.50 of A < 0.6, so no
#  intersection: 0, union 200, 0, 0.
###########################################################
echo -e "    jaccard.new.t012...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n0\t200\t0\t0\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -f 0.6 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t013: -f keeps an overlap of exactly the fraction
#  Same input as t012 with -f 0.5: 50 150 0.333333 1.
###########################################################
echo -e "    jaccard.new.t013...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n50\t150\t0.333333\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -f 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t014: -F tests the fraction of B
#  A chr1:0-100, B chr1:50-250 (200 bp): 50 bp = 0.25 of B < 0.3:
#  intersection 0, union 100 + 200 = 300.
###########################################################
echo -e "    jaccard.new.t014...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t250\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n0\t300\t0\t0\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -F 0.3 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t015: -f with -r requires the fraction of both
#  Same input as t014 with -f 0.5 -r: 0.50 of A passes, 0.25 of B
#  fails: 0 300 0 0.
###########################################################
echo -e "    jaccard.new.t015...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t250\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n0\t300\t0\t0\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -f 0.5 -r > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t016: -e accepts either fraction
#  Same input as t014 with -f 0.9 -F 0.2 -e: 0.25 of B passes.
#  Intersection 50, union 300 - 50 = 250, jaccard 0.2, 1 pair.
###########################################################
echo -e "    jaccard.new.t016...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t250\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n50\t250\t0.2\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -f 0.9 -F 0.2 -e > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t017: -f is tested against the merged interval
#  A chr1:0-50 and chr1:40-100 merge to chr1:0-100. B chr1:30-90
#  covers 60 bp = 0.60 of the merged A, so -f 0.6 keeps it
#  (against chr1:0-50 alone it would be 20/50 = 0.40).
#  Intersection 60, union 100 + 60 - 60 = 100, jaccard 0.6.
###########################################################
echo -e "    jaccard.new.t017...\c"
printf 'chr1\t0\t50\nchr1\t40\t100\n' > a.tmp.bed
printf 'chr1\t30\t90\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n60\t100\t0.6\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -f 0.6 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t018: -s merges and intersects per strand
#  A chr1:0-100 (+). B chr1:50-150 (-) and chr1:60-80 (+) are kept
#  apart because their strands differ. Only the + B intersects A:
#  20 bp. Union 100 + (100 + 20) - 20 = 200, jaccard 0.1, 1 pair.
#  pins the current output: with -s each file is merged per strand;
#  the docs say only that hits must be on the same strand.
###########################################################
echo -e "    jaccard.new.t018...\c"
printf 'chr1\t0\t100\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t50\t150\tb1\t0\t-\nchr1\t60\t80\tb2\t0\t+\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n20\t200\t0.1\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -s > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t019: -s drops unknown-strand records from every
#  column
#  A and B each have chr1:0-100 (+) and chr1:200-300 (.). Without
#  -s this is 200 200 1 2; with -s only the + records count:
#  100 100 1 1. pins the current output: with -s unknown-strand
#  records are removed from both files, union included; the docs
#  say only that hits must be on the same strand.
###########################################################
echo -e "    jaccard.new.t019...\c"
printf 'chr1\t0\t100\ta1\t0\t+\nchr1\t200\t300\ta2\t0\t.\n' > a.tmp.bed
printf 'chr1\t0\t100\tb1\t0\t+\nchr1\t200\t300\tb2\t0\t.\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n100\t100\t1\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -s > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t023: -split drops a B that falls between blocks
#  A chr1:0-100 has blocks 0-10 and 90-100 (20 bp); B chr1:40-60
#  lies in the gap. Intersection 0; union uses the block length
#  of A: 20 + 20 = 40.
###########################################################
echo -e "    jaccard.new.t023...\c"
printf 'chr1\t0\t100\ta1\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,\n' > a.tmp.bed
printf 'chr1\t40\t60\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n0\t40\t0\t0\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -split > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t024: -split counts only the bases inside blocks
#  Same A as t023; B chr1:5-15 overlaps block 0-10 by 5 bp.
#  Intersection 5, union 20 + 10 - 5 = 25, jaccard 0.2, 1 pair.
###########################################################
echo -e "    jaccard.new.t024...\c"
printf 'chr1\t0\t100\ta1\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,\n' > a.tmp.bed
printf 'chr1\t5\t15\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n5\t25\t0.2\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -split > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t025: jaccard has no -N option
#  -N belongs to closest; jaccard rejects it as unrecognized.
###########################################################
echo -e "    jaccard.new.t025...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\n' > b.tmp.bed
printf '***** ERROR: Unrecognized parameter: -N *****\nexit 1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -N > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t026: -header prints the query's header before the
#  table
#  Same intervals as t003.
###########################################################
echo -e "    jaccard.new.t026...\c"
printf '# header\nchr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\n' > b.tmp.bed
printf '# header\nintersection\tunion\tjaccard\tn_intersections\n50\t150\t0.333333\t1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  jaccard.new.t027: a GFF query uses 0-based coordinates
#  GFF chr1 1-100 is chr1:0-100; with B chr1:50-150 this is t003:
#  50 150 0.333333 1.
###########################################################
echo -e "    jaccard.new.t027...\c"
printf 'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=a1\n' > a.tmp.gff
printf 'chr1\t50\t150\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n50\t150\t0.333333\t1\n' > exp
$BT jaccard -a a.tmp.gff -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  jaccard.new.t028: adjacent VCF records merge
#  SNVs at POS 1 and 2 are chr1:0-1 and chr1:1-2, book-ended, so
#  they merge to chr1:0-2. B chr1:0-2: intersection 2, union
#  2 + 2 - 2 = 2, jaccard 1, 1 pair. Unmerged records would give 2
#  pairs; reading POS as 0-based (chr1:1-3) would give 1 3 0.333333 1.
###########################################################
echo -e "    jaccard.new.t028...\c"
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\trs1\tA\tG\t50\tPASS\t.\nchr1\t2\trs2\tC\tT\t50\tPASS\t.\n' > a.tmp.vcf
printf 'chr1\t0\t2\n' > b.tmp.bed
printf 'intersection\tunion\tjaccard\tn_intersections\n2\t2\t1\t1\n' > exp
$BT jaccard -a a.tmp.vcf -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  jaccard.new.t029: an unsorted query is an error
###########################################################
echo -e "    jaccard.new.t029...\c"
printf 'chr1\t100\t200\nchr1\t0\t50\n' > a.tmp.bed
printf 'chr1\t50\t150\n' > b.tmp.bed
printf 'Error: Sorted input specified, but the file a.tmp.bed has the following out of order record\nchr1\t0\t50\nexit 1\n' > exp
$BT jaccard -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

[[ $FAILURES -eq 0 ]] || exit 1;
