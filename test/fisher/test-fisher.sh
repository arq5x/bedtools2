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

echo -e "    fisher.t1...\c"
echo \
"# Number of query intervals: 3
# Number of db intervals: 2
# Number of overlaps: 2
# Number of possible intervals (estimated): 34
# phyper(2 - 1, 3, 34 - 3, 2, lower.tail=F)
# Contingency Table Of Counts
#_________________________________________
#           |  in -b       | not in -b    |
#     in -a | 2            | 1            |
# not in -a | 0            | 31           |
#_________________________________________
# p-values for fisher's exact test
left	right	two-tail	ratio
1	0.0053476	0.0053476	inf" > exp
$BT fisher -a a.bed -b b.bed -g t.500.genome > obs
check obs exp
rm obs exp


echo -e "    fisher.t2...\c"
echo \
"# Number of query intervals: 3
# Number of db intervals: 2
# Number of overlaps: 2
# Number of possible intervals (estimated): 4
# phyper(2 - 1, 3, 4 - 3, 2, lower.tail=F)
# Contingency Table Of Counts
#_________________________________________
#           |  in -b       | not in -b    |
#     in -a | 2            | 1            |
# not in -a | 0            | 1            |
#_________________________________________
# p-values for fisher's exact test
left	right	two-tail	ratio
1	0.5	1	inf" > exp
$BT fisher -a a.bed -b b.bed -g t.60.genome > obs
check obs exp
rm obs exp


echo -e "    fisher.t3...\c"
echo \
"# Number of query intervals: 4
# Number of db intervals: 2
# Number of overlaps: 3
# Number of possible intervals (estimated): 4
# phyper(3 - 1, 4, 4 - 4, 2, lower.tail=F)
# Contingency Table Of Counts
#_________________________________________
#           |  in -b       | not in -b    |
#     in -a | 3            | 1            |
# not in -a | 0            | 0            |
#_________________________________________
# p-values for fisher's exact test
left	right	two-tail	ratio
1	1	1	-nan" > exp
$BT fisher -a a_merge.bed -b b.bed -g t.60.genome > obs
check obs exp
rm obs exp


echo -e "    fisher.t4...\c"
echo \
"# Number of query intervals: 3
# Number of db intervals: 2
# Number of overlaps: 2
# Number of possible intervals (estimated): 4
# phyper(2 - 1, 3, 4 - 3, 2, lower.tail=F)
# Contingency Table Of Counts
#_________________________________________
#           |  in -b       | not in -b    |
#     in -a | 2            | 1            |
# not in -a | 0            | 1            |
#_________________________________________
# p-values for fisher's exact test
left	right	two-tail	ratio
1	0.5	1	inf" > exp
$BT fisher -a a_merge.bed -b b.bed -g t.60.genome -m > obs
check obs exp
rm obs exp
[[ $FAILURES -eq 0 ]] || exit 1;

echo -e "    fisher.t5...\c"
$BT fisher -b test.bed -a tumor.gff -g dm6.fai > exp
DIR="${TMPDIR:-/tmp}/Users/mvandenb/src/genomic_features_bardin_lab/build/modencode/dm6/bed/"
LONG_PATH="$DIR/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAaS-\(phospho-S2\)_.bed"
mkdir -p "$DIR"
cp test.bed "$LONG_PATH"
$BT fisher -b "$LONG_PATH" -a tumor.gff -g dm6.fai > obs
check obs exp
rm obs exp
[[ $FAILURES -eq 0 ]] || exit 1;


# Regression test for issue 954
echo -e "    fisher.t6...\c"
echo \
"# Number of query intervals: 5
# Number of db intervals: 9
# Number of overlaps: 2
# Number of possible intervals (estimated): 74
# phyper(2 - 1, 5, 74 - 5, 9, lower.tail=F)
# Contingency Table Of Counts
#_________________________________________
#           |  in -b       | not in -b    |
#     in -a | 2            | 3            |
# not in -a | 7            | 62           |
#_________________________________________
# p-values for fisher's exact test
left	right	two-tail	ratio
0.98864	0.10898	0.10898	5.905" > exp
$BT fisher -a issue954_a.bed -b issue954_b.bed -g issue954.genome > obs
check obs exp
rm obs exp

###########################################################
#
# NEW TESTS: table, estimate, p-values and options
#
# fisher counts A records (qc), B records (dc) and overlapping
# A/B pairs (n11); n12 = max(0, qc-n11), n21 = max(0, dc-n11).
# The number of possible intervals is
# max(n11+n12+n21, floor(genome / ((1 + A bp/qc) + (1 + B bp/dc))))
# and n22 is that minus n11+n12+n21.
# p-values: hypergeometric P(k) = C(K,k)C(N-K,n-k)/C(N,n) with
# N = n11+n12+n21+n22, K = n11+n12, n = n11+n21; left = P(k <= n11),
# right = P(k >= n11), two-tail = sum of P(k) <= P(n11). Printed %.5g;
# ratio (n11/n12)/(n21/n22) printed %.3f.
# Derivation of each p-value, standard library only (set a, b, c, d
# to n11, n12, n21, n22 and run the lines below without "#  "):
#  python3 - <<'EOF'
#  from math import comb as C
#  a, b, c, d = 2, 1, 1, 41
#  N = a + b + c + d; K = a + b; n = a + c
#  P = {k: C(K, k) * C(N - K, n - k) / C(N, n)
#       for k in range(max(0, K + n - N), min(K, n) + 1)}
#  print("%.5g %.5g %.5g" % (sum(p for k, p in P.items() if k <= a),
#        sum(p for k, p in P.items() if k >= a),
#        min(1, sum(p for p in P.values() if p <= P[a] * (1 + 1e-7)))))
#  EOF
#
###########################################################

###########################################################
#  fisher.new.t001: the table, estimate and p-values for a basic case
#  Genome 1000 bp. A: 3 records of 10 bp (union 30); B: 3 of 10 bp
#  (union 30). Overlapping pairs: A1-B1 and A2-B2, so n11 = 2.
#  Mean sizes: 1 + 30/3 = 11 each, sum 22; possible intervals
#  floor(1000/22) = 45. n12 = 3-2 = 1, n21 = 3-2 = 1,
#  n22 = 45-2-1-1 = 41. N = 45, K = 3, n = 3, C(45,3) = 14190:
#  P(0) = 11480, P(1) = 2583, P(2) = 126, P(3) = 1 (/14190).
#  left = 14189/14190 = 0.99993; right = two-tail = 127/14190 =
#  0.00895; ratio (2/1)/(1/41) = 82.000.
#  pins the current output: %.5g p-values, %.3f ratio.
###########################################################
echo -e "    fisher.new.t001...\c"
printf 'chr1\t1000\n' > g.tmp
printf 'chr1\t0\t10\nchr1\t100\t110\nchr1\t500\t510\n' > a.tmp.bed
printf 'chr1\t5\t15\nchr1\t105\t115\nchr1\t800\t810\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 3' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 2' \
	'# Number of possible intervals (estimated): 45' \
	'# phyper(2 - 1, 3, 45 - 3, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 2            | 1            |' \
	'# not in -a | 1            | 41           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.99993\t0.00895\t0.00895\t82.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t002: without -m, overlapping records count separately
#  A: chr1:0-10, 5-15, 500-510, 600-610 (4 records, 40 bp).
#  B: chr1:5-15, 12-20, 800-810, 900-910 (4 records, 38 bp).
#  Pairs: 0-10/5-15, 5-15/5-15, 5-15/12-20, so n11 = 3, n12 = 1,
#  n21 = 1. Mean sizes 1 + 40/4 = 11 and 1 + 38/4 = 10.5, sum 21.5;
#  floor(1000/21.5) = 46, n22 = 46-5 = 41. N = 46, K = 4, n = 4:
#  right = two-tail = (C(4,3)C(42,1) + 1)/C(46,4) = 169/163185 =
#  0.0010356; left = 1 - 1/163185 = 0.99999; ratio 3/(1/41) = 123.
###########################################################
echo -e "    fisher.new.t002...\c"
printf 'chr1\t1000\n' > g.tmp
printf 'chr1\t0\t10\nchr1\t5\t15\nchr1\t500\t510\nchr1\t600\t610\n' > a.tmp.bed
printf 'chr1\t5\t15\nchr1\t12\t20\nchr1\t800\t810\nchr1\t900\t910\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 4' \
	'# Number of db intervals: 4' \
	'# Number of overlaps: 3' \
	'# Number of possible intervals (estimated): 46' \
	'# phyper(3 - 1, 4, 46 - 4, 4, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 3            | 1            |' \
	'# not in -a | 1            | 41           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.99999\t0.0010356\t0.0010356\t123.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t003: -m merges overlapping records in both files first
#  Same input as t002. Merged A: chr1:0-15, 500-510, 600-610 (3, 35
#  bp); merged B: chr1:5-20, 800-810, 900-910 (3, 35 bp). One pair:
#  n11 = 1, n12 = 2, n21 = 2. Mean sizes 1 + 35/3 each, sum 25.33;
#  floor(1000/25.33) = 39, n22 = 39-5 = 34. N = 39, K = 3, n = 3,
#  C(39,3) = 9139: P(0) = C(36,3) = 7140, so right = 1999/9139 =
#  0.21873, left = 1 - C(3,2)C(36,1)/9139 - 1/9139 = 9030/9139 =
#  0.98807; ratio (1/2)/(2/34) = 8.500. -m is given last: its handler
#  skips the argument that follows it.
###########################################################
echo -e "    fisher.new.t003...\c"
printf 'chr1\t1000\n' > g.tmp
printf 'chr1\t0\t10\nchr1\t5\t15\nchr1\t500\t510\nchr1\t600\t610\n' > a.tmp.bed
printf 'chr1\t5\t15\nchr1\t12\t20\nchr1\t800\t810\nchr1\t900\t910\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 3' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 1' \
	'# Number of possible intervals (estimated): 39' \
	'# phyper(1 - 1, 3, 39 - 3, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 1            | 2            |' \
	'# not in -a | 2            | 34           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.98807\t0.21873\t0.21873\t8.500' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp -m > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t004: -s counts only same-strand overlaps
#  Genome 2000 bp. A: a1 (+), a2 (+), a3 (-); B: b1 (+), b2 (-),
#  b3 (-), each 10 bp, each a overlapping its b. Same strand: a1-b1
#  and a3-b3, n11 = 2, n12 = 1, n21 = 1. Mean sizes sum 22,
#  floor(2000/22) = 90, n22 = 86. N = 90, K = 3, n = 3:
#  right = two-tail = (C(3,2)C(87,1) + 1)/C(90,3) = 262/117480 =
#  0.0022302; left = 1 - 1/117480 = 0.99999; ratio 2/(1/86) = 172.
###########################################################
echo -e "    fisher.new.t004...\c"
printf 'chr1\t2000\n' > g.tmp
printf 'chr1\t0\t10\ta1\t0\t+\nchr1\t100\t110\ta2\t0\t+\nchr1\t500\t510\ta3\t0\t-\n' > a.tmp.bed
printf 'chr1\t5\t15\tb1\t0\t+\nchr1\t105\t115\tb2\t0\t-\nchr1\t505\t515\tb3\t0\t-\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 3' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 2' \
	'# Number of possible intervals (estimated): 90' \
	'# phyper(2 - 1, 3, 90 - 3, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 2            | 1            |' \
	'# not in -a | 1            | 86           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.99999\t0.0022302\t0.0022302\t172.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp -s > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t005: one A overlapping two B records counts two
#  overlaps
#  Genome 10000 bp. A: chr1:0-100, 200-300, 800-900, 950-960 (4, 310
#  bp); B: chr1:90-110, 95-105, 250-450, 600-680 (4, 310 bp). Mean
#  sizes 1 + 310/4 = 78.5 each, sum 157; floor(10000/157) = 63.
#  Overlapping pairs (fraction of A / of B):
#  0-100/90-110 10 bp (0.10/0.50), 0-100/95-105 5 bp (0.05/0.50),
#  200-300/250-450 50 bp (0.50/0.25).
#  n11 = 3, n12 = 1, n21 = 1, n22 = 63-5 = 58. N = 63, K = 4, n = 4:
#  right = two-tail = (C(4,3)C(59,1) + 1)/C(63,4) = 237/595665 =
#  0.00039787; left = 1 - 1/595665 = 1 at 5 digits; ratio 3/(1/58)
#  = 174.
###########################################################
echo -e "    fisher.new.t005...\c"
printf 'chr1\t10000\n' > g.tmp
printf 'chr1\t0\t100\nchr1\t200\t300\nchr1\t800\t900\nchr1\t950\t960\n' > a.tmp.bed
printf 'chr1\t90\t110\nchr1\t95\t105\nchr1\t250\t450\nchr1\t600\t680\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 4' \
	'# Number of db intervals: 4' \
	'# Number of overlaps: 3' \
	'# Number of possible intervals (estimated): 63' \
	'# phyper(3 - 1, 4, 63 - 4, 4, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 3            | 1            |' \
	'# not in -a | 1            | 58           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'1\t0.00039787\t0.00039787\t174.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t006: -f counts pairs covering the fraction of A
#  Same input as t005. With -f 0.5 only 200-300/250-450 (0.50 of A)
#  counts: n11 = 1, n12 = 3, n21 = 3, n22 = 63-7 = 56. N = 63, K = 4,
#  n = 4: right = 1 - C(59,4)/C(63,4) = 1 - 455126/595665 = 0.23594,
#  left = 0.98237 (1 - P(k >= 2)); ratio (1/3)/(3/56) = 6.222.
###########################################################
echo -e "    fisher.new.t006...\c"
printf 'chr1\t10000\n' > g.tmp
printf 'chr1\t0\t100\nchr1\t200\t300\nchr1\t800\t900\nchr1\t950\t960\n' > a.tmp.bed
printf 'chr1\t90\t110\nchr1\t95\t105\nchr1\t250\t450\nchr1\t600\t680\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 4' \
	'# Number of db intervals: 4' \
	'# Number of overlaps: 1' \
	'# Number of possible intervals (estimated): 63' \
	'# phyper(1 - 1, 4, 63 - 4, 4, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 1            | 3            |' \
	'# not in -a | 3            | 56           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.98237\t0.23594\t0.23594\t6.222' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp -f 0.5 > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t007: -F counts pairs covering the fraction of B
#  Same input as t005. With -F 0.5 the two pairs with 0.50 of B
#  count: n11 = 2, n12 = 2, n21 = 2, n22 = 63-6 = 57. N = 63, K = 4,
#  n = 4: right = P(k >= 2) = 0.017632, left = 0.9996; ratio
#  (2/2)/(2/57) = 28.500.
###########################################################
echo -e "    fisher.new.t007...\c"
printf 'chr1\t10000\n' > g.tmp
printf 'chr1\t0\t100\nchr1\t200\t300\nchr1\t800\t900\nchr1\t950\t960\n' > a.tmp.bed
printf 'chr1\t90\t110\nchr1\t95\t105\nchr1\t250\t450\nchr1\t600\t680\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 4' \
	'# Number of db intervals: 4' \
	'# Number of overlaps: 2' \
	'# Number of possible intervals (estimated): 63' \
	'# phyper(2 - 1, 4, 63 - 4, 4, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 2            | 2            |' \
	'# not in -a | 2            | 57           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.9996\t0.017632\t0.017632\t28.500' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp -F 0.5 > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t008: -f with -r needs the fraction of A and of B
#  Same input as t005. With -f 0.5 -r no pair has 0.50 of both:
#  n11 = 0, n12 = 4, n21 = 4, n22 = 63-8 = 55. N = 63, K = 4, n = 4:
#  left = P(0) = C(59,4)/C(63,4) = 455126/595665 = 0.76406, right = 1,
#  two-tail = 1 (P(0) is the most likely table); ratio 0/(4/55) =
#  0.000.
###########################################################
echo -e "    fisher.new.t008...\c"
printf 'chr1\t10000\n' > g.tmp
printf 'chr1\t0\t100\nchr1\t200\t300\nchr1\t800\t900\nchr1\t950\t960\n' > a.tmp.bed
printf 'chr1\t90\t110\nchr1\t95\t105\nchr1\t250\t450\nchr1\t600\t680\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 4' \
	'# Number of db intervals: 4' \
	'# Number of overlaps: 0' \
	'# Number of possible intervals (estimated): 63' \
	'# phyper(0 - 1, 4, 63 - 4, 4, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 0            | 4            |' \
	'# not in -a | 4            | 55           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.76406\t1\t1\t0.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp -f 0.5 -r > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t009: -e counts pairs passing either fraction
#  Same input as t005. With -f 0.5 -F 0.5 -e all three pairs pass one
#  of the two, so the table is that of t005.
###########################################################
echo -e "    fisher.new.t009...\c"
printf 'chr1\t10000\n' > g.tmp
printf 'chr1\t0\t100\nchr1\t200\t300\nchr1\t800\t900\nchr1\t950\t960\n' > a.tmp.bed
printf 'chr1\t90\t110\nchr1\t95\t105\nchr1\t250\t450\nchr1\t600\t680\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 4' \
	'# Number of db intervals: 4' \
	'# Number of overlaps: 3' \
	'# Number of possible intervals (estimated): 63' \
	'# phyper(3 - 1, 4, 63 - 4, 4, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 3            | 1            |' \
	'# not in -a | 1            | 58           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'1\t0.00039787\t0.00039787\t174.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp -f 0.5 -F 0.5 -e > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t010: -split counts block overlaps and block lengths
#  Genome 1000 bp. A1 chr1:0-100 has blocks 0-10 and 90-100 (20 bp);
#  A2 chr1:300-310 has one block (10 bp): A union 30. B: chr1:40-60
#  (in A1's gap), 95-105 (in A1's second block), 500-510: 3 records,
#  40 bp. Only 95-105 counts: n11 = 1, n12 = 1, n21 = 2. Mean sizes
#  1 + 30/2 = 16 and 1 + 40/3, sum 30.33; floor(1000/30.33) = 32,
#  n22 = 32-4 = 28. N = 32, K = 2, n = 3: P(0) = C(30,3)/C(32,3) =
#  4060/4960, P(2) = 30/4960; right = 900/4960 = 0.18145, left =
#  1 - 30/4960 = 0.99395; ratio (1/1)/(2/28) = 14.000.
###########################################################
echo -e "    fisher.new.t010...\c"
printf 'chr1\t1000\n' > g.tmp
printf 'chr1\t0\t100\ta1\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,\nchr1\t300\t310\ta2\t0\t+\t300\t310\t0\t1\t10,\t0,\n' > a.tmp.bed
printf 'chr1\t40\t60\nchr1\t95\t105\nchr1\t500\t510\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 2' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 1' \
	'# Number of possible intervals (estimated): 32' \
	'# phyper(1 - 1, 2, 32 - 2, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 1            | 1            |' \
	'# not in -a | 2            | 28           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.99395\t0.18145\t0.18145\t14.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp -split > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t011: -exclude removes the merged excluded length from
#  the genome size
#  Same A and B as t001; the excluded chr1:850-900 and 880-1000 merge
#  to 150 bp, so the genome is 1000-150 = 850 bp. floor(850/22) = 38,
#  n22 = 38-4 = 34. N = 38, K = 3, n = 3: right = two-tail =
#  (C(3,2)C(35,1) + 1)/C(38,3) = 106/8436 = 0.012565, left = 1 -
#  1/8436 = 0.99988; ratio 2/(1/34) = 68.000. -exclude is accepted
#  by the parser but missing from the usage text.
#  The expectation is correct, but master computes it only while
#  BedFile::_total_flattened_length starts at zero: no constructor
#  initializes it (bedFile.h), so this test fails when fresh heap
#  memory is not zero (e.g. MALLOC_PERTURB_).
###########################################################
echo -e "    fisher.new.t011...\c"
printf 'chr1\t1000\n' > g.tmp
printf 'chr1\t0\t10\nchr1\t100\t110\nchr1\t500\t510\n' > a.tmp.bed
printf 'chr1\t5\t15\nchr1\t105\t115\nchr1\t800\t810\n' > b.tmp.bed
printf 'chr1\t850\t900\nchr1\t880\t1000\n' > ex.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 3' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 2' \
	'# Number of possible intervals (estimated): 38' \
	'# phyper(2 - 1, 3, 38 - 3, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 2            | 1            |' \
	'# not in -a | 1            | 34           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.99988\t0.012565\t0.012565\t68.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -exclude ex.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp ex.tmp.bed a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t012: the genome size sums every chromosome, and B
#  records on a chromosome A lacks are counted
#  Genome chr1 500 + chr2 500 + chr3 200 = 1200 bp. A: chr1:0-10,
#  chr1:100-110, chr2:0-10; B: chr1:5-15, chr2:5-15, chr3:0-10 (3 db
#  intervals). Pairs: 2. floor(1200/22) = 54, n22 = 54-4 = 50.
#  N = 54, K = 3, n = 3: right = two-tail = (3*51 + 1)/C(54,3) =
#  154/24804 = 0.0062087, left = 1 - 1/24804 = 0.99996; ratio
#  2/(1/50) = 100.000.
###########################################################
echo -e "    fisher.new.t012...\c"
printf 'chr1\t500\nchr2\t500\nchr3\t200\n' > g.tmp
printf 'chr1\t0\t10\nchr1\t100\t110\nchr2\t0\t10\n' > a.tmp.bed
printf 'chr1\t5\t15\nchr2\t5\t15\nchr3\t0\t10\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 3' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 2' \
	'# Number of possible intervals (estimated): 54' \
	'# phyper(2 - 1, 3, 54 - 3, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 2            | 1            |' \
	'# not in -a | 1            | 50           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.99996\t0.0062087\t0.0062087\t100.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t013: n12 is clamped at 0 when one A has more overlaps
#  than there are A records
#  A: chr1:0-100 (1 record); B: chr1:10-20, 30-40, 500-510 (30 bp).
#  Pairs: 2, so n12 = max(0, 1-2) = 0, n21 = 1. Mean sizes 101 and
#  11, sum 112; floor(1000/112) = 8 >= 3, n22 = 8-3 = 5. N = 8,
#  K = 2, n = 3, C(8,3) = 56: P(0) = 20, P(1) = 30, P(2) = 6 (/56);
#  left = 1, right = two-tail = 6/56 = 0.10714; ratio (2/0)/(1/5) is
#  inf.
###########################################################
echo -e "    fisher.new.t013...\c"
printf 'chr1\t1000\n' > g.tmp
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\nchr1\t30\t40\nchr1\t500\t510\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 1' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 2' \
	'# Number of possible intervals (estimated): 8' \
	'# phyper(2 - 1, 1, 8 - 1, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 2            | 0            |' \
	'# not in -a | 1            | 5            |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'1\t0.10714\t0.10714\tinf' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t014: the estimate never falls below n11+n12+n21
#  Genome 40 bp. A: 3 records of 5 bp; B: 3 of 5 bp; 2 pairs. Mean
#  sizes 6 each, sum 12; floor(40/12) = 3 < 2+1+1 = 4, so the
#  estimate is 4 and n22 = 0. N = 4, K = 3, n = 3: P(2) = 3/4,
#  P(3) = 1/4; left = 0.75, right = 1, two-tail = 1; ratio
#  (2/1)/(1/0) = 0.000.
###########################################################
echo -e "    fisher.new.t014...\c"
printf 'chr1\t40\n' > g.tmp
printf 'chr1\t0\t5\nchr1\t10\t15\nchr1\t20\t25\n' > a.tmp.bed
printf 'chr1\t2\t7\nchr1\t12\t17\nchr1\t30\t35\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 3' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 2' \
	'# Number of possible intervals (estimated): 4' \
	'# phyper(2 - 1, 3, 4 - 3, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 2            | 1            |' \
	'# not in -a | 1            | 0            |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.75\t1\t1\t0.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t015: book-ended records do not overlap; a depleted
#  table gives a small left p-value
#  Genome 100 bp. Five 1 bp A records, each book-ended by a 1 bp B
#  record: no pairs. Mean sizes 2 each, sum 4; floor(100/4) = 25.
#  n11 = 0, n12 = 5, n21 = 5, n22 = 15. N = 25, K = 5, n = 5,
#  C(25,5) = 53130: P(0..5) = 15504, 24225, 11400, 1900, 100, 1.
#  left = 15504/53130 = 0.29181, right = 1, two-tail = (15504 +
#  11400 + 1900 + 100 + 1)/53130 = 0.54404; ratio 0.000.
###########################################################
echo -e "    fisher.new.t015...\c"
printf 'chr1\t100\n' > g.tmp
printf 'chr1\t0\t1\nchr1\t10\t11\nchr1\t20\t21\nchr1\t30\t31\nchr1\t40\t41\n' > a.tmp.bed
printf 'chr1\t1\t2\nchr1\t11\t12\nchr1\t21\t22\nchr1\t31\t32\nchr1\t41\t42\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 5' \
	'# Number of db intervals: 5' \
	'# Number of overlaps: 0' \
	'# Number of possible intervals (estimated): 25' \
	'# phyper(0 - 1, 5, 25 - 5, 5, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 0            | 5            |' \
	'# not in -a | 5            | 15           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.29181\t1\t0.54404\t0.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t016: GFF records are converted to 0-based before sizes
#  and overlaps are measured
#  GFF 1-10, 101-110, 501-510 are chr1:0-10, 100-110, 500-510, the A
#  of t001, so the output is that of t001. Reading the GFF coordinates
#  as 0-based (9 bp records) would change the estimate.
###########################################################
echo -e "    fisher.new.t016...\c"
printf 'chr1\t1000\n' > g.tmp
printf 'chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=a1\nchr1\tsrc\tgene\t101\t110\t.\t+\t.\tID=a2\nchr1\tsrc\tgene\t501\t510\t.\t+\t.\tID=a3\n' > a.tmp.gff
printf 'chr1\t5\t15\nchr1\t105\t115\nchr1\t800\t810\n' > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 3' \
	'# Number of db intervals: 3' \
	'# Number of overlaps: 2' \
	'# Number of possible intervals (estimated): 45' \
	'# phyper(2 - 1, 3, 45 - 3, 3, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 2            | 1            |' \
	'# not in -a | 1            | 41           |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'0.99993\t0.00895\t0.00895\t82.000' > exp
$BT fisher -a a.tmp.gff -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.gff b.tmp.bed

###########################################################
#  fisher.new.t017: a very small p-value is printed in exponent form
#  Genome 1000000 bp. Ten 10 bp A records and ten 10 bp B records;
#  the first nine pairs overlap. Mean sizes sum 22,
#  floor(1000000/22) = 45454; n11 = 9, n12 = 1, n21 = 1,
#  n22 = 45454-11 = 45443. N = 45454, K = 10, n = 10: right =
#  two-tail = (C(10,9)C(45444,1) + 1)/C(45454,10) = 4.3848e-35,
#  left = 1; ratio 9/(1/45443) = 408987.000.
#  pins the current output: %.5g switches to exponent form.
###########################################################
echo -e "    fisher.new.t017...\c"
printf 'chr1\t1000000\n' > g.tmp
printf 'chr1\t%d\t%d\n' 0 10 1000 1010 2000 2010 3000 3010 4000 4010 5000 5010 6000 6010 7000 7010 8000 8010 9000 9010 > a.tmp.bed
printf 'chr1\t%d\t%d\n' 5 15 1005 1015 2005 2015 3005 3015 4005 4015 5005 5015 6005 6015 7005 7015 8005 8015 9500 9510 > b.tmp.bed
printf '%b\n' \
	'# Number of query intervals: 10' \
	'# Number of db intervals: 10' \
	'# Number of overlaps: 9' \
	'# Number of possible intervals (estimated): 45454' \
	'# phyper(9 - 1, 10, 45454 - 10, 10, lower.tail=F)' \
	'# Contingency Table Of Counts' \
	'#_________________________________________' \
	'#           |  in -b       | not in -b    |' \
	'#     in -a | 9            | 1            |' \
	'# not in -a | 1            | 45443        |' \
	'#_________________________________________' \
	'# p-values for fisher'\''s exact test' \
	'left\tright\ttwo-tail\tratio' \
	'1\t4.3848e-35\t4.3848e-35\t408987.000' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  fisher.new.t018: -g is required
###########################################################
echo -e "    fisher.new.t018...\c"
printf 'chr1\t0\t10\nchr1\t100\t110\nchr1\t500\t510\n' > a.tmp.bed
printf 'chr1\t5\t15\nchr1\t105\t115\nchr1\t800\t810\n' > b.tmp.bed
printf 'ERROR*****: specify -g genome file*****\n\nexit 1\n' > exp
$BT fisher -a a.tmp.bed -b b.tmp.bed > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

[[ $FAILURES -eq 0 ]] || exit 1;
