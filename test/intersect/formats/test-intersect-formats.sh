set -e;
echo -e \
"\n###########################################################
#
#  INTERSECT FORMATS
#
###########################################################\n"

BT=${BT-../../../bin/bedtools}

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


[[ $FAILURES -eq 0 ]] || exit 1;
