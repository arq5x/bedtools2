set -e;
###########################################################
#
#  Unit tests for sampleFile program
#
############################################################

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
#
# Create a 1,000 record BED6 file named "mainFile.bed".
# Give it a one line header, to test the -header option.
#
###########################################################
echo "#This is the mainFile from which samples will be taken." > mainFile.bed
$BT random -l 1000 -n 1000 -g human.hg19.genome >> mainFile.bed



###########################################################
#  Test that help is printed when no args are given
############################################################
echo -e "    sample.t01...\c"
echo \
"***** ERROR: No input file given. Exiting. *****" > exp
$BT sample 2>&1 > /dev/null | tail -1 > obs
check obs exp
rm obs


###########################################################
#  Test that we throw an error for unrecognized arguments
############################################################
echo -e "    sample.new.t02...\c"
echo "***** ERROR: Unrecognized parameter: -wrongArg *****" > exp
$BT sample -i mainFile.bed -wrongArg 2>&1 > /dev/null |  tail -1 > obs
check obs exp
rm obs exp

###########################################################
#  Test that we throw an error when no input file was given.
############################################################
echo -e "    sample.new.t03...\c"
echo "***** ERROR: No input file given. Exiting. *****" > exp;
$BT sample -n 10 2>&1 > /dev/null | tail -1 > obs
check obs exp
rm obs exp

###########################################################
#  Test that we throw an error for -i without input file
############################################################
echo -e "    sample.new.t04...\c"
echo "***** ERROR: -i option given, but no input file specified. *****" > exp
$BT sample -i 2>&1 > /dev/null | tail -1 > obs
check obs exp
rm obs exp


###########################################################
#  Test that we throw an error for -n given without 
#  number of output records sepcified.
############################################################
echo -e "    sample.new.t05...\c"
echo "***** ERROR: -n option given, but no number of output records specified. *****" > exp
$BT sample -n 2>&1 > /dev/null | tail -1 > obs
check obs exp
rm obs exp


###########################################################
#  Test that we throw an error when num output records
#  exceeds records in file.
############################################################
echo -e "    sample.new.t06...\c"
echo "***** ERROR: Input file has fewer records than the requested number of output records. *****" > exp
$BT sample -i mainFile.bed 2>&1 > /dev/null | head -2 | tail -1 > obs
check obs exp
rm obs exp


###########################################################
#  Test that we get the requested number of records
############################################################
echo -e "    sample.new.t07...\c"
echo "10" > exp
# BSD wc adds leading whitespace, and BSD sed -i param requires extension
$BT sample -i mainFile.bed -n 10 | wc -l | sed 's/^[ \t]*//' > obs
check obs exp
rm obs exp


###########################################################
#  Test that the -seed option gives consistent results
############################################################
echo -e "    sample.new.t08...\c"
$BT sample -i mainFile.bed -n 50 -seed 4 > obs
$BT sample -i mainFile.bed -n 50 -seed 4 > exp
check obs exp
rm obs exp


###########################################################
#  Test that -header option gives header
############################################################
echo -e "    sample.new.t09...\c"
echo "#This is the mainFile from which samples will be taken." > exp
$BT sample -i mainFile.bed -n 10 -header | head -1 > obs
check obs exp
rm obs exp




rm mainFile.bed





[[ $FAILURES -eq 0 ]] || exit 1;
