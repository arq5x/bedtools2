BT=${BT-../../bin/bedtools}

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

###########################################################
#  Test a basic file
###########################################################
echo "    spacing.t01...\c"
echo \
"chr1	20	30	.
chr1	25	40	0
chr1	50	80	10
chr1	75	100	0
chr1	105	110	5
chr2	115	130	.
chr2	120	160	0
chr2	170	180	10" > exp
$BT spacing -i a.bed  > obs
check obs exp
rm obs exp

