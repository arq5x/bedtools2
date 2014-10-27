BT=${BT-../../bin/bedtools}

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

samtools view -Sb test.sam > test.bam 2> /dev/null

$BT bamtofastq -i test.bam -fq test.fq -fq2 test.fq2 2> /dev/null

check test.fq golden.fq
check test.fq2 golden.fq2

rm test.bam test.fq test.fq2

