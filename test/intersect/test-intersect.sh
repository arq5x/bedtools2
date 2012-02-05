check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

echo "    intersect.t1 self intersect...\c"
../../bin/bedtools intersect -a a.bed -b a.bed > obs
check obs intersect.t1.exp
