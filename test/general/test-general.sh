BT=../../bin/bedtools

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

###########################################################
#  Enforce non-negative coordinates
###########################################################
echo "    general.t01...\c"
echo \
"chr1	1	10
chr1	-1	10" | $BT merge -i - 2> obs
echo "Error: malformed BED entry at line 2. Start Coordinate detected that is < 0. Exiting." > exp
check obs exp
rm obs exp

###########################################################
#  Enforce start <= end
###########################################################
echo "    general.t02...\c"
echo \
"chr1	1	2
chr1	10	5" | $BT merge -i - 2> obs
echo "Error: malformed BED entry at line 2. Start was greater than end. Exiting." > exp
check obs exp
rm obs exp

###########################################################
#  Enforce integer coordinates
###########################################################
echo "    general.t03...\c"
echo \
"chr1	.	2" | $BT merge -i - 2> obs
echo "Unexpected file format.  Please use tab-delimited BED, GFF, or VCF. Perhaps you have non-integer starts or ends at line 1?" > exp
check obs exp
rm obs exp


###########################################################
#  Enforce integer coordinates
###########################################################
echo "    general.t04...\c"
echo \
"chr1	.	2" | $BT merge -i - 2> obs
echo "Unexpected file format.  Please use tab-delimited BED, GFF, or VCF. Perhaps you have non-integer starts or ends at line 1?" > exp
check obs exp
rm obs exp


###########################################################
#  Enforce tab-separated files
###########################################################
echo "    general.t05...\c"
echo \
"chr1 1 2" | $BT merge -i - 2> obs
echo "It looks as though you have less than 3 columns at line: 1.  Are you sure your files are tab-delimited?" > exp
check obs exp
rm obs exp

