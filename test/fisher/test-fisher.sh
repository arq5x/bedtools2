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

echo "    fisher.t1...\c"
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


echo "    fisher.t2...\c"
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


echo "    fisher.t3...\c"
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


echo "    fisher.t4...\c"
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
