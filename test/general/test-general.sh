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


###########################################################
#  Fail on non-existent files.
###########################################################
echo "    general.t06...\c"
$BT merge -i idontexist.bed 2> obs
echo "Error: The requested file (idontexist.bed) could not be opened. Error message: (No such file or directory). Exiting!" > exp
check obs exp
rm obs exp


###########################################################
#  Don't fail on existent, yet empty files.
###########################################################
echo "    general.t07...\c"
$BT merge -i empty.bed > obs
touch exp
check obs exp
rm obs exp


###########################################################
#  Process gzipped files.
###########################################################
echo "    general.t08...\c"
$BT merge -i non-empty.bed.gz > obs
echo "chr1	10	21" > exp
check obs exp
rm obs exp


