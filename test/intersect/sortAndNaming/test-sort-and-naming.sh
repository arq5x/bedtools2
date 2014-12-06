echo -e \
"\n###########################################################
#  
#  CHROMOSOME SORT ORDER AND NAMING CONVENTIONS 
#
###########################################################\n"


BT=${BT-../../../bin/bedtools}

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

###########################################################
#  Test query against DB with reverse order
############################################################
echo "    intersect.t01...\c"
echo \
"ERROR: chromomsome sort ordering for file sdb1.bed is inconsistent with other files. Record was:
chr2	10	20" > exp
$BT intersect -a sq1.bed -b sdb1.bed -sorted 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs


###########################################################
#  Test query with 2dbs, one of which is out of order
############################################################
echo "    intersect.t02...\c"
echo \
"ERROR: chromomsome sort ordering for file db2_numBackwards.bed is inconsistent with other files. Record was:
chr3	60	100" > exp
$BT intersect -a q1a_num.bed -b db1_num.bed db2_numBackwards.bed -sorted 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs


###########################################################
#  Test query with 3dbs, one of which is out of order
############################################################
echo "    intersect.t03...\c"
echo \
"ERROR: chromomsome sort ordering for file db3_numBackwards.bed is inconsistent with other files. Record was:
chr1	10	20" > exp
$BT intersect -a q1a_num.bed -b db1_num.bed db2_num.bed db3_numBackwards.bed -sorted 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs


###########################################################
#  Test query with 2 dbs, one of which has a chrom 
#  that the query does not
############################################################
echo "    intersect.t04...\c"
echo \
"ERROR: Database file db1_num.bed contains chromosome chr3, but the query file does not.
       Please re-reun with the -g option for a genome file.
       See documentation for details." >exp
$BT intersect -a q1_num.bed -b db1_num.bed db2_num.bed -sorted 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs

###########################################################
#  Test query with db that is missing "chr" in one of it's
#  chroms.
############################################################
echo "    intersect.t05...\c"
echo \
"ERROR: File db1_noChr.bed has inconsistent naming convention for record:
2	50	80" > exp
$BT intersect -a q1_num.bed -b db1_noChr.bed  -sorted 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs

###########################################################
#  Test query with db that has leading zero in it's chrom 
#  names.
############################################################
echo "    intersect.t06...\c"
echo \
"ERROR: File db1_leadingZero.txt has a record where naming convention (leading zero) is inconsistent with other files:
chr01	10	20" > exp
$BT intersect -a q1_num.bed -b db1_leadingZero.txt -sorted 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs

###########################################################
#  Test that leading zeroes are allowed if they appear after
# an underscore
############################################################
echo "    intersect.t07...\c"
echo \
"chr1	10	20
chr1	80	100
chr2	50	80
chr2	100	120
chr10	5	50
chr10	80	120
chr11	20	60
chr11	80	120
chr12	10	50
chr12	60	90
chr1_gl0003	20	80" > exp
$BT intersect -a q1_gls.bed -b q1_gls.bed  -sorted > obs
check exp obs
rm exp obs
