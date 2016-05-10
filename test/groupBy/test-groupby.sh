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
#  Test that -n option is shown as deperecated
###########################################################
#echo "    merge.t2...\c"
#echo "***** ERROR: -n option is deprecated. Please see the documentation for the -c and -o column operation options. *****" > exp
#$BT merge -i a.bed -n 2>&1 > /dev/null | head -2 | tail -1 > obs
#check obs exp
#rm obs exp


###########################################################
#  Test basic grouping
###########################################################
echo "    groupby.t1...\c"
echo \
"chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	45
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	8" > exp
$BT groupby -i values3.header.bed -c 5  > obs
check obs exp
rm obs exp

###########################################################
#  Test case insensitive grouping works
###########################################################
echo "    groupby.t2...\c"
echo \
"chr1	0	10	10
cHr1	10	20	5
Chr1	11	21	5
chR1	20	30	45
Chr1	120	130	1
CHr3	0	10	1
cHR3	10	20	2
CHR3	20	30	3
chr3	120	130	8" > exp
$BT groupby -i values3_case.header.bed -c 5 -ignorecase > obs
check obs exp
rm obs exp


###########################################################
#  Test -full option (print all columns, not just grouped
#  ones)
###########################################################
echo "    groupby.t3...\c"
echo \
"chr1	0	10	a1	10	+	10
chr1	10	20	a2	5	+	5
chr1	11	21	a3	5	+	5
chr1	20	30	a4	15	+	45
chr1	120	130	a7	1	+	1
chr3	0	10	a8	1	+	1
chr3	10	20	a9	2	+	2
chr3	20	30	a10	3	+	3
chr3	120	130	a11	4	+	8" > exp
$BT groupby -i values3.header.bed -c 5 -full > obs
check obs exp
rm obs exp



###########################################################
#  Test -inheader option
###########################################################
echo "    groupby.t4...\c"
echo \
"chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	45
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	8" > exp
$BT groupby -i values3.header.bed -c 5 -inheader > obs
check obs exp
rm obs exp

###########################################################
#  Test -inheader option when header not marked by
#  recognized char
###########################################################
echo "    groupby.t5...\c"
echo \
"chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	45
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	8" > exp
$BT groupby -i values3.unmarked_header.bed -c 5 -inheader > obs
check obs exp
rm obs exp

###########################################################
#  Test -inheader option when no header present will skip
# first line
###########################################################
echo "    groupby.t6...\c"
echo \
"chr1	10	20	5
chr1	11	21	5
chr1	20	30	45
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	8" > exp
$BT groupby -i values3.no_header.bed -c 5 -inheader > obs
check obs exp
rm obs exp

###########################################################
#  Test -outheader option will work automatically, even
# without -inheader, if header has normally marked start char.
###########################################################
echo "    groupby.t7...\c"
echo \
"#chrom	start	end	A	B	C
chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	45
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	8" > exp
$BT groupby -i values3.header.bed -c 5 -outheader > obs
check obs exp
rm obs exp

###########################################################
#  Test that unmarked header will be included by default.
###########################################################
echo "    groupby.t8...\c"
echo \
"chrom	start	end	B
chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	15
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	4" > exp
$BT groupby -i values3.unmarked_header.bed -c 5 -o distinct > obs
check obs exp
rm obs exp

###########################################################
#  Test that -outheader does nothing with unmarked header
###########################################################
echo "    groupby.t9...\c"
echo \
"col_1	col_2	col_3	col_4	col_5	col_6
chrom	start	end	B
chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	15
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	4" > exp
$BT groupby -i values3.unmarked_header.bed -c 5 -o distinct -outheader > obs
check obs exp
rm obs exp

###########################################################
#  Test that -header works with unmarked header
###########################################################
echo "    groupby.t10...\c"
echo \
"chrom	start	end	A	B	C
chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	15
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	4" > exp
$BT groupby -i values3.unmarked_header.bed -c 5 -o distinct -header > obs
check obs exp
rm obs exp

###########################################################
#  Test that -header works normally with normal header
###########################################################
echo "    groupby.t11...\c"
echo \
"#chrom	start	end	A	B	C
chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	45
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	8" > exp
$BT groupby -i values3.header.bed -c 5  -header > obs
check obs exp
rm obs exp


###########################################################
#  Test a BedPlus file (7 fields)
###########################################################
echo "    groupby.t12...\c"
echo \
"chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	45
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	8" > exp
$BT groupby -i values3.7fields.header.bed -c 5  > obs
check obs exp
rm obs exp


###########################################################
#  Test noPosPlus file (8 fields, not starting with 
# chr, starte, end
###########################################################
echo "    groupby.t13...\c"
echo \
"chr1	0	10	10
chr1	10	20	5
chr1	11	21	5
chr1	20	30	45
chr1	120	130	1
chr3	0	10	1
chr3	10	20	2
chr3	20	30	3
chr3	120	130	8" > exp
$BT groupby -g 2-4 -i noPosvalues.header.bed -c 6 > obs
check obs exp
rm obs exp

###########################################################
#  Test noPosPlus file with mof columns (iterated and range)
###########################################################
echo "    groupby.t14...\c"
echo \
"0	10	chr1	10
10	20	chr1	5
11	21	chr1	5
20	30	chr1	45
120	130	chr1	1
0	10	chr3	1
10	20	chr3	2
20	30	chr3	3
120	130	chr3	8" > exp
$BT groupby -g 3-4,2 -i noPosvalues.header.bed -c 6 > obs
check obs exp
rm obs exp

###########################################################
#  Test that only the groupBy tool may use 
# non-positional records
###########################################################
echo "    groupby.t15...\c"
echo \
"ERROR: file noPosvalues.header.bed has non positional records, which are only valid for the groupBy tool." > exp
$BT merge  -i noPosvalues.header.bed 2>&1 >/dev/null | cat - > obs
check obs exp
rm obs exp

###########################################################
#  Test a VCF file
###########################################################
echo "    groupby.t16...\c"
echo \
"19	G	70.9
19	C	33.71
19	A	21.2" > exp
$BT groupby -i a_vcfSVtest.vcf -g 1,4 -o mean -c 6 > obs
check obs exp
rm obs exp


###########################################################
#  Test a BAM file
###########################################################
echo "    groupby.t17...\c"
echo \
"None	chr2L	118.75" > exp
$BT groupby -i gdc.bam -g 1,3 -c 4 -o mean > obs
check obs exp
rm obs exp

###########################################################
#  Test a single column of data
###########################################################
echo "    groupby.t18...\c"
echo \
"chr1	chr1,chr1,chr1" > exp
cut -f 1 test.bed | $BT groupby -g 1 -i - -c 1 -o collapse > obs
check obs exp
rm obs exp
