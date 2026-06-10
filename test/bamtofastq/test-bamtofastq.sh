set -e;
BT=${BT-../../bin/bedtools}
htsutil=${htsutil-../htsutil}

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

$htsutil samtobam test2.sam test2.bam

$BT bamtofastq -i test2.bam -fq test2.fq -fq2 test2.fq2 2> /dev/null

check test2.fq golden.fq
check test2.fq2 golden.fq2

rm test2.bam test2.fq test2.fq2

###########################################################
# Test: single-end with secondary/supplementary alignments
# should only output primary alignments (no duplicates)
###########################################################
echo "    bamtofastq.t2...\c"
$htsutil samtobam secondary_single.sam secondary_single.bam
$BT bamtofastq -i secondary_single.bam -fq secondary_single_out.fq 2> /dev/null
check secondary_single_out.fq golden_secondary_single.fq
rm secondary_single.bam secondary_single_out.fq

###########################################################
# Test: paired-end with secondary alignments
# should only output primary alignment pairs (no duplicates)
###########################################################
echo "    bamtofastq.t3...\c"
$htsutil samtobam secondary_paired.sam secondary_paired.bam
$BT bamtofastq -i secondary_paired.bam -fq secondary_paired_out.fq -fq2 secondary_paired_out.fq2 2> /dev/null
check secondary_paired_out.fq golden_secondary_paired.fq
check secondary_paired_out.fq2 golden_secondary_paired.fq2
rm secondary_paired.bam secondary_paired_out.fq secondary_paired_out.fq2

[[ $FAILURES -eq 0 ]] || exit 1;
