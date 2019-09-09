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

samtools view -Sb test2.sam > test2.bam 2> /dev/null

$BT bamtofastq -i test2.bam -fq test2.fq -fq2 test2.fq2 2> /dev/null

check test2.fq golden.fq
check test2.fq2 golden.fq2

rm test2.bam test2.fq test2.fq2

[[ $FAILURES -eq 0 ]] || exit 1;
