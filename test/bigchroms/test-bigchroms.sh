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
############################################################
echo -e "    bigchroms.t01...\c"
$BT intersect -sorted -a abig.bed -b abig.bed > obs
check obs abig.bed
rm obs
