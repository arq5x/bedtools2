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
[[ $FAILURES -eq 0 ]] || exit 1;
