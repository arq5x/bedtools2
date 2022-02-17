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

# cat a.bed
# chr1	10	20
# chr1	30	40
# chr1	40	50
# chr1	45	100

###########################################################
#  Test a basic merge; one interval should be un-merged, 
#  the other two should be merged.
###########################################################
echo -e "    merge.t1...\c"
echo \
"chr1	10	20
chr1	30	100" > exp
$BT merge -i a.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test that -n option is shown as deperecated
###########################################################
echo -e "    merge.t2...\c"
echo "***** ERROR: -n option is deprecated. Please see the documentation for the -c and -o column operation options. *****" > exp
$BT merge -i a.bed -n 2>&1 > /dev/null | tail -1 > obs
check obs exp
rm obs exp


###########################################################
#  Test the counting of merged intervals. (-n)
###########################################################
echo -e "    merge.t3...\c"
echo \
"chr1	10	20	1
chr1	30	100	3" > exp
$BT merge -i a.bed -c 1 -o count > obs
check obs exp
rm obs exp


###########################################################
#  Test that -nms option is deprecated
###########################################################
echo -e "    merge.t4...\c"
echo "***** ERROR: -nms option is deprecated. Please see the documentation for the -c and -o column operation options. *****" > exp
$BT merge -i a.bed -nms 2>&1 > /dev/null  | tail -1 > obs
check obs exp
rm obs exp

###########################################################
#  Test the listing of names from merged intervals.
###########################################################
echo -e "    merge.t5...\c"
echo \
"chr1	10	20	a1
chr1	30	100	a2,a3,a4" > exp
$BT merge -i a.names.bed -c 4 -o collapse > obs
check obs exp
rm obs exp

###########################################################
# collapsed list of the names, and sum of the scores
###########################################################
echo -e "    merge.t6...\c"
echo \
"chr1	10	20	a1	1
chr1	30	100	a2,a3,a4	9
chr2	10	20	a1	5
chr2	30	40	a2	6
chr2	42	100	a3,a4	15" > exp
$BT merge -i a.full.bed -c 4,5  -o collapse,sum > obs
check obs exp
rm obs exp

###########################################################
# count intervals and sum of scores
###########################################################
echo -e "    merge.t7...\c"
echo \
"chr1	10	20	1	1
chr1	30	100	3	9
chr2	10	20	1	5
chr2	30	40	1	6
chr2	42	100	2	15" > exp
$BT merge -i a.full.bed -c 5 -o count,sum> obs
check obs exp
rm obs exp

###########################################################
# count, collapsed names, and sum of scores
###########################################################
echo -e "    merge.t8...\c"
echo \
"chr1	10	20	a1	1	1
chr1	30	100	a2,a3,a4	9	3
chr2	10	20	a1	5	1
chr2	30	40	a2	6	1
chr2	42	100	a3,a4	15	2" > exp
$BT merge -i a.full.bed -c 4,5,4 -o collapse,sum,count > obs
check obs exp
rm obs exp

###########################################################
# stranded merge, show sign, collapsed names, sum of
# scores, and count
###########################################################
echo -e "    merge.t9a...\c"
echo \
"chr1	10	20	a1	1	1
chr1	30	40	a2	2	1
chr1	40	50	a3	3	1
chr1	45	100	a4	4	1
chr2	10	20	a1	5	1
chr2	30	40	a2	6	1
chr2	42	50	a3	7	1
chr2	45	100	a4	8	1" > exp
$BT merge -i a.full.bed -s -c 4,5,6 -o collapse,sum,count > obs
check obs exp
rm obs exp

###########################################################
# stranded merge, show sign, collapsed names, sum of
# scores, and count
###########################################################
echo -e "    merge.t9b...\c"
echo \
"chr1	10	20	a1	1	+
chr1	30	40	a2	2	+
chr1	40	50	a3	3	-
chr1	45	100	a4	4	+
chr2	10	20	a1	5	+
chr2	30	40	a2	6	+
chr2	42	50	a3	7	+
chr2	45	100	a4	8	-" > exp
$BT merge -i a.full.bed -s -c 4,5,6 -o collapse,sum,collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test the use of a custom delimiter for -delim option
###########################################################
echo -e "    merge.t10...\c"
echo \
"chr1	10	20	a1
chr1	30	100	a2|a3|a4" > exp
$BT merge -i a.names.bed -delim "|" -c 4 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test that stranded merge not allowed with VCF
###########################################################
VCF_FILE=testA.vcf
echo -e "    merge.t11...\c"
echo "***** ERROR: stranded merge not supported for VCF file $VCF_FILE. *****" >exp
$BT merge -i $VCF_FILE -s 2>&1 > /dev/null | tail -1 > obs
check exp obs
rm obs exp

###########################################################
#  Test that column ops not allowed with BAM if col greater
#  than 11.
#  
# EDIT: This test has been moved to test #35, after the 
# other bam column tests. 
###########################################################
echo -e "    merge.t12...\c"
#echo "***** ERROR: Requested column 12, but database file fullFields.bam only has fields 1 - 11." > exp
#$BT merge -i fullFields.bam -c 12 -o sum 2>&1 > /dev/null | head -3 | tail -1 > obs
#check exp obs
#rm obs exp
echo ok

###########################################################
#  Test that VCF input gives BED3 output
###########################################################
echo -e "    merge.t13...\c"
echo \
"chr1	30859	30860
chr1	69269	69270
chr1	69510	69511
chr1	874815	874816
chr1	879675	879676
chr1	935491	935492
chr1	1334051	1334057
chr1	31896607	31896608" > exp
$BT merge -i testA.vcf > obs
check exp obs
rm obs exp

###########################################################
#  Test that GFF input gives BED3 output
###########################################################
echo -e "    merge.t14...\c"
echo \
"chr22	9999999	10001000
chr22	10009999	10010100
chr22	10019999	10025000" > exp
$BT merge -i a.gff > obs
check exp obs
rm obs exp

###########################################################
#  Test that stranded merge where some records have
#  unknown strand works correctly
###########################################################
echo -e "    merge.t15...\c"
echo \
"chr1	10	80
chr1	20	90
chr2	20	60
chr2	25	80" > exp
$BT merge -i mixedStrands.bed -s  > obs
check exp obs
rm obs exp

###########################################################
#  Test that stranded merge with unknown records works
#  correctly, forward strand only
###########################################################
echo -e "    merge.t16...\c"
echo \
"chr1	10	80
chr2	20	60" > exp
$BT merge -i mixedStrands.bed -S + > obs
check exp obs
rm obs exp

###########################################################
#  Test that stranded merge with unknown records works
#  correctly, reverse strand only
###########################################################
echo -e "    merge.t17...\c"
echo \
"chr1	20	90
chr2	25	80" > exp
$BT merge -i mixedStrands.bed -S - > obs
check exp obs
rm obs exp

###########################################################
#  Test that merge with specified strand does not allowed
#  other characters besides + or -.
###########################################################
echo -e "    merge.t18...\c"
echo "***** ERROR: -S option must be followed by + or -. *****" > exp
$BT merge -i mixedStrands.bed -S . -c 6 -o distinct 2>&1 > /dev/null | tail -1 >obs
check exp obs
rm obs exp


###########################################################
#  Test that sort order is enforced
###########################################################
echo -e "    merge.t19...\c"
echo \
"Error: Sorted input specified, but the file unsorted.bed has the following out of order record
chr1	9	30	2" > exp
$BT merge -i unsorted.bed 2>&1 > /dev/null | tail -2  >obs
check exp obs
rm obs exp

###########################################################
#  Test that chrom change is handled correctly
###########################################################
echo -e "    merge.t20...\c"
echo \
"chr1	9	30
chr1	100	110
chr2	11	20" > exp
$BT merge -i b.bed > obs
check exp obs
rm exp obs

###########################################################
#  Test that a merged BAM file only gives BED3 output
###########################################################
echo -e "    merge.t21...\c"
echo \
"chr1	10	20
chr1	30	100
chr2	10	20
chr2	30	40
chr2	42	100" > exp
$BT merge -i a.full.bed > obs
check exp obs
rm exp obs

###########################################################
#  Test that precision is correct
###########################################################
echo -e "    merge.t22...\c"
echo \
"chr2L	1	54	0.05
chr2L	65	128	0.33
chr2L	129	180	0.04
chr2L	193	317	-0.125
chr2L	321	375	-0.07
chr2L	385	448	-0.11
chr2L	449	502	0.4
chr2L	513	570	0.48
chr2L	577	635	-0.24" > exp
$BT merge -i precisionTest.bed -c 5 -o mean > obs
check obs exp
rm obs exp

###########################################################
#  Test that numeric ops on non-numeric columns
#  are allowed, but produce a warning and null
#  value result. 
###########################################################
echo -e "    merge.t23a...\c"
echo \
"chr1	10	20	.
chr1	30	100	." > expOut
$BT merge -i a.names.bed -c 4 -o sum 2>&1 > obsOut | cat - > obsErr
check obsOut expOut
rm expOut obsOut


###########################################################
#  Just check that the warning message from the previous
#  test was correct.
###########################################################
echo -e "    merge.t23b...\c"
echo \
" ***** WARNING: Non numeric value a1 in 4.
 ***** WARNING: Non numeric value a4 in 4." > expErr
check obsErr expErr
rm obsErr expErr


###########################################################
#
#  Test that we can get fields from a BAM file
#
###########################################################

###########################################################
#  Test bam column 1
###########################################################
echo -e "    merge.t24...\c"
$BT merge -i fullFields.bam -c 1 -o collapse > obs
check obs bamCol1Collapse.txt
rm obs

###########################################################
#  Test bam column 2 gives an error
###########################################################
echo -e "    merge.t25...\c"
echo \
"***** ERROR: Requested column 2 of a BAM file, which is the Flags field." > exp
$BT merge -i fullFields.bam -c 2 -o collapse 2>&1 > /dev/null | head -3 | tail -1 > obs
check exp obs
rm obs exp

###########################################################
#  Test bam column 3
###########################################################
echo -e "    merge.t26...\c"
$BT merge -i fullFields.bam -c 3 -o collapse > obs
check obs bamCol3Collapse.txt
rm obs

###########################################################
#  Test bam column 4
###########################################################
echo -e "    merge.t27...\c"
$BT merge -i fullFields.bam -c 4 -o mean > obs
check obs bamCol4Mean.txt
rm obs



###########################################################
#  Test bam column 5
###########################################################
echo -e "    merge.t28...\c"
$BT merge -i fullFields.bam -c 5 -o mean > obs
check obs bamCol5Mean.txt
rm obs

###########################################################
#  Test bam column 6
###########################################################
echo -e "    merge.t29...\c"
$BT merge -i fullFields.bam -c 6 -o collapse > obs
check obs bamCol6Collapse.txt
rm obs

###########################################################
#  Test bam column 7
###########################################################
echo -e "    merge.t30...\c"
$BT merge -i fullFields.bam -c 7 -o collapse > obs
check obs bamCol7Collapse.txt
rm obs

###########################################################
#  Test bam column 8
###########################################################
echo -e "    merge.t31...\c"
$BT merge -i fullFields.bam -c 8 -o mean > obs
check obs bamCol8Mean.txt
rm obs

###########################################################
#  Test bam column 9
###########################################################
echo -e "    merge.t32...\c"
$BT merge -i fullFields.bam -c 9 -o mean > obs
check obs bamCol9Mean.txt
rm obs

###########################################################
#  Test bam column 10
###########################################################
echo -e "    merge.t33...\c"
$BT merge -i fullFields.bam -c 10 -o collapse > obs
check obs bamCol10Collapse.txt
rm obs

###########################################################
#  Test bam column 11
###########################################################
echo -e "    merge.t34...\c"
$BT merge -i fullFields.bam -c 11 -o collapse > obs
check obs bamCol11Collapse.txt
rm obs

###########################################################
#  Test that column ops not allowed with BAM if col greater
#  than 11.
###########################################################
echo -e "    merge.t35...\c"
echo "***** ERROR: Requested column 12, but database file fullFields.bam only has fields 1 - 11." > exp
$BT merge -i fullFields.bam -c 12 -o sum 2>&1 > /dev/null | head -3 | tail -1 > obs
check exp obs
rm obs exp

###########################################################
#  Test col ops behavior on bam file for missing values,
#  i.e, getting the mate reference when there is no mate.
#  Be sure null value is printed
###########################################################
echo -e "    merge.t36...\c"
echo \
"chr1	10	20	.
chr1	30	100	.,.,.
chr2	10	20	.
chr2	30	40	.
chr2	42	100	.,." >exp
$BT merge -i a.full.bam  -c 7 -o collapse > obs
check exp obs
rm obs exp




###########################################################
#
#  Test new -iobuf option
#  
###########################################################


###########################################################
#  Test -iobuf expects an argument
###########################################################
echo -e "    merge.t37...\c"
echo "***** ERROR: -iobuf option given, but size of input buffer not specified. *****" >exp
$BT merge -i a.bed -iobuf 2>&1 > /dev/null |  tail -1 > obs
check obs exp
rm obs exp

###########################################################
#  Test -iobuf allows only suffixes K/M/G
###########################################################
echo -e "    merge.t38...\c"
echo \
"***** ERROR: Unrecognized memory buffer size suffix 'L' given. *****" > exp
$BT merge -i a.bed -iobuf 20L 2>&1 > /dev/null | tail -1 > obs
check obs exp
rm obs exp

###########################################################
#  Test -iobuf doesn't allow a buffer size below 8 bytes.
###########################################################
echo -e "    merge.t39...\c"
echo \
"***** ERROR: specified buffer size is too small. *****" > exp
$BT merge -i a.bed -iobuf 7 2>&1 > /dev/null | tail -1 > obs
check exp obs
rm exp obs

###########################################################
#  Test -iobuf doesn't allow non-numeric arguments
###########################################################
echo -e "    merge.t40...\c"
echo \
"***** ERROR: argument passed to -iobuf is not numeric. *****" > exp
$BT merge -i a.bed -iobuf beerM 2>&1 > /dev/null | tail -1 > obs
check exp obs
rm exp obs

###########################################################
#  Test -iobuf allows correct argument with suffix
###########################################################
echo -e "    merge.t41...\c"
echo \
"chr1	10	20
chr1	30	100" > exp
$BT merge -i a.bed -iobuf 128M > obs
check exp obs
rm exp obs

###########################################################
#  Test -iobuf allows correct argument without suffix
###########################################################
echo -e "    merge.t42...\c"
echo \
"chr1	10	20
chr1	30	100" > exp
$BT merge -i a.bed -iobuf 8192 > obs
check exp obs
rm exp obs

###########################################################
#  Test that scientific notation is allowed for coordinates
###########################################################
echo -e "    merge.t43...\c"
echo \
"chr1	800	830" > exp
$BT merge -i expFormat.bed > obs
check exp obs
rm obs exp


###########################################################
#  Test that struct vars in VCF get correct length
###########################################################
echo -e "    merge.t44a...\c"
echo \
"19	252805	257416
19	260364	261044
19	265133	265691
19	265985	266386" > exp
$BT merge -i vcfSVtest.vcf > obs
check exp obs
rm obs exp

###########################################################
#  Test that struct vars in VCF get correct length
###########################################################
echo -e "    merge.t44b...\c"
if [[ -f vcfSVtest.2.vcf ]]; then
    echo \
        "19	252805	297416" > exp
    $BT merge -i vcfSVtest.2.vcf > obs
    check exp obs
    rm obs exp
else
    echo "skipped - could not find vcfSVtest.2.vcf";
fi

###########################################################
#  Test that stdin is used by default
###########################################################
echo -e "    merge.t45...\c"
echo \
"chr1	10	20
chr1	30	100" >exp
cat a.bed | $BT merge > obs
check exp obs
rm obs exp


###########################################################
#  Test that precision default is high enough for 
#  formatting not to give scientific notation
###########################################################
echo -e "    merge.t46...\c"
echo \
"chr1	5333587	5344172	5344172
chr1	5481008	5484749	16454247
chr1	6763278	6766882	6766882" > exp
$BT merge -i precisionTest2.bed -c 8 -o sum> obs
check exp obs
rm obs exp


###########################################################
#  Test that user can specify a lower precision
###########################################################
echo -e "    merge.t47...\c"
echo \
"chr1	5333587	5344172	5.3442e+06
chr1	5481008	5484749	1.6454e+07
chr1	6763278	6766882	6.7669e+06" > exp
$BT merge -i precisionTest2.bed -c 8 -o sum -prec 5 > obs
check exp obs
rm obs exp

###########################################################
#  Test stranded merge with bedPlus files that have strand
###########################################################
echo -e "    merge.t48...\c"
echo \
"chr1	10000	25000" > exp
$BT merge -i bug254_d.bed -s -d 200 > obs
check exp obs
rm obs exp

###########################################################
#  Test stranded merge with bedPlus files that have strand
###########################################################
echo -e "    merge.t49...\c"
echo \
"chr1	10000	20000
chr1	20100	25000" > exp
$BT merge -i bug254_e.bed -s -d 200 > obs
check exp obs
rm obs exp


###########################################################
#  Test the chained GZIP file - Regression test for bug #975
###########################################################
echo -e "    merge.t50...\c"
echo \
"1	10000	20000
5	55554	66666"	>	exp
$BT merge -i chained.bed.gz > obs
check exp obs
rm obs exp

[[ $FAILURES -eq 0 ]] || exit 1;
