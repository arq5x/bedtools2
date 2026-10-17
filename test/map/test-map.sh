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
#  Test defaults
############################################################
echo -e "    map.t01...\c"
echo \
"chr1	0	100	30
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	6
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test sum
############################################################
echo -e "    map.t02...\c"
echo \
"chr1	0	100	30
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	6
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed -o sum > obs
check obs exp
rm obs exp


###########################################################
#  Test count
############################################################
echo -e "    map.t03...\c"
echo \
"chr1	0	100	3
chr1	100	200	1
chr2	0	100	0
chr2	100	200	0
chr3	0	100	3
chr3	100	200	1" > exp
$BT map -a ivls.bed -b values.bed -o count > obs
check obs exp
rm obs exp


###########################################################
#  Test mean
############################################################
echo -e "    map.t04...\c"
echo \
"chr1	0	100	10
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	2
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed -o mean > obs
check obs exp
rm obs exp

###########################################################
#  Test max
############################################################
echo -e "    map.t05...\c"
echo \
"chr1	0	100	15
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	3
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed -o max > obs
check obs exp
rm obs exp

###########################################################
#  Test min
############################################################
echo -e "    map.t06...\c"
echo \
"chr1	0	100	5
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	1
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values.bed -o min > obs
check obs exp
rm obs exp

###########################################################
#  Test mode
############################################################
echo -e "    map.t07...\c"
echo \
"chr1	0	100	5
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	1
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values2.bed -o mode > obs
check obs exp
rm obs exp

###########################################################
#  Test anti-mode
############################################################
echo -e "    map.t08...\c"
echo \
"chr1	0	100	10
chr1	100	200	1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	1
chr3	100	200	4" > exp
$BT map -a ivls.bed -b values2.bed -o antimode > obs
check obs exp
rm obs exp


###########################################################
#  Test column extraction from BEDPLUS
############################################################
echo -e "    map.t09...\c"
echo \
"chr1	0	100	1,2,3,4,5,-6
chr1	100	200	7
chr2	0	100	.
chr2	100	200	.
chr3	0	100	8,9,-10
chr3	100	200	11,12" > exp
$BT map -a ivls.bed -b values4.bed -c 7 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test min
############################################################
echo -e "    map.t10...\c"
echo \
"chr1	0	100	-6
chr1	100	200	7
chr2	0	100	.
chr2	100	200	.
chr3	0	100	-10
chr3	100	200	11" > exp
$BT map -a ivls.bed -b values4.bed -c 7 -o min > obs
check obs exp
rm obs exp

###########################################################
#  Test absmin
############################################################
echo -e "    map.t11...\c"
echo \
"chr1	0	100	1
chr1	100	200	7
chr2	0	100	.
chr2	100	200	.
chr3	0	100	8
chr3	100	200	11" > exp
$BT map -a ivls.bed -b values4.bed -c 7 -o absmin > obs
check obs exp
rm obs exp

###########################################################
#  Test max
############################################################
echo -e "    map.t12...\c"
echo \
"chr1	0	100	5
chr1	100	200	7
chr2	0	100	.
chr2	100	200	.
chr3	0	100	9
chr3	100	200	12" > exp
$BT map -a ivls.bed -b values4.bed -c 7 -o max > obs
check obs exp
rm obs exp

###########################################################
#  Test absmax
############################################################
echo -e "    map.t13...\c"
echo \
"chr1	0	100	6
chr1	100	200	7
chr2	0	100	.
chr2	100	200	.
chr3	0	100	10
chr3	100	200	12" > exp
$BT map -a ivls.bed -b values4.bed -c 7 -o absmax > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t14...\c"
echo \
"chr1	0	100	chr1,chr1,chr1,chr1,chr1
chr1	100	200	chr1
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 1 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t15...\c"
echo \
"chr1	0	100	hg19_ccdsGene,hg19_ccdsGene,hg19_ccdsGene,hg19_ccdsGene,hg19_ccdsGene
chr1	100	200	hg19_ccdsGene
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 2 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t16...\c"
echo \
"chr1	0	100	start_codon,CDS,exon,CDS,exon
chr1	100	200	exon
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 3 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t17...\c"
echo \
"chr1	0	100	1,2,8,9,40
chr1	100	200	40
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 4 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t18...\c"
echo \
"chr1	0	100	9,11,20,17,200
chr1	100	200	200
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 5 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t19...\c"
echo \
"chr1	0	100	0.000000,0.000000,0.000000,0.000000,0.000000
chr1	100	200	0.000000
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 6 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t20...\c"
echo \
"chr1	0	100	+,+,+,+,+
chr1	100	200	+
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 7 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t21...\c"
echo \
"chr1	0	100	.,0,.,2,.
chr1	100	200	.
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 8 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test GFF column extraction
############################################################
echo -e "    map.t22..\c"
echo \
"chr1	0	100	gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";,gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";,gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";,gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";,gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";
chr1	100	200	gene_id \"CCDS30744.1\"; transcript_id \"CCDS30744.1\";
chr2	0	100	.
chr2	100	200	.
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.gff -c 9 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t23..\c"
echo \
"chr1	0	100	chr1,chr1,chr1
chr1	100	200	chr1,chr1
chr2	0	100	.
chr2	100	200	chr2,chr2
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 1 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t24..\c"
echo \
"chr1	0	100	10,15,20
chr1	100	200	110,130
chr2	0	100	.
chr2	100	200	110,130
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 2 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t25..\c"
echo \
"chr1	0	100	rs6054257,.,rs6040355
chr1	100	200	.,microsat1
chr2	0	100	.
chr2	100	200	.,microsat1
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 3 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t26..\c"
echo \
"chr1	0	100	G,T,A
chr1	100	200	T,GTC
chr2	0	100	.
chr2	100	200	T,GTC
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 4 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t27..\c"
echo \
"chr1	0	100	A,A,G,T
chr1	100	200	.,G,GTCT
chr2	0	100	.
chr2	100	200	.,G,GTCT
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 5 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t28..\c"
echo \
"chr1	0	100	29,3,67
chr1	100	200	47,50
chr2	0	100	.
chr2	100	200	47,50
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 6 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t29..\c"
echo \
"chr1	0	100	PASS,q10,PASS
chr1	100	200	PASS,PASS
chr2	0	100	.
chr2	100	200	PASS,PASS
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 7 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t30..\c"
echo \
"chr1	0	100	NS=3;DP=14;AF=0.5;DB;H2,NS=3;DP=11;AF=0.017,NS=2;DP=10;AF=0.333,0.667;AA=T;DB
chr1	100	200	NS=3;DP=13;AA=T,NS=3;DP=9;AA=G
chr2	0	100	.
chr2	100	200	NS=3;DP=13;AA=T,NS=3;DP=9;AA=G
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 8 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t31..\c"
echo \
"chr1	0	100	GT:GQ:DP:HQ,GT:GQ:DP:HQ,GT:GQ:DP:HQ
chr1	100	200	GT:GQ:DP:HQ,GT:GQ:DP
chr2	0	100	.
chr2	100	200	GT:GQ:DP:HQ,GT:GQ:DP
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 9 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t32..\c"
echo \
"chr1	0	100	0|0:48:1:51,51,0|0:49:3:58,50,1|2:21:6:23,27
chr1	100	200	0|0:54:7:56,60,0/1:35:4
chr2	0	100	.
chr2	100	200	0|0:54:7:56,60,0/1:35:4
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -c 10 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test VCF column extraction
############################################################
echo -e "    map.t33..\c"
echo \
"
*****
***** ERROR: Requested column 15, but database file test.vcf only has fields 1 - 12." > exp
$BT map -a ivls.bed -b test.vcf -c 15 -o collapse 2>&1 > /dev/null | head -3> obs
check obs exp
rm obs exp

###########################################################
#  Test -header
############################################################
echo -e "    map.t33..\c"
echo \
"#header
chr1	0	100	0|0:48:1:51,51,0|0:49:3:58,50,1|2:21:6:23,27
chr1	100	200	0|0:54:7:56,60,0/1:35:4
chr2	0	100	.
chr2	100	200	0|0:54:7:56,60,0/1:35:4
chr3	0	100	.
chr3	100	200	." > exp
$BT map -a ivls.bed -b test.vcf -header -c 10 -o collapse > obs
check obs exp
rm obs exp

###########################################################
#  Test -null
############################################################
echo -e "    map.t33..\c"
echo \
"chr1	0	100	0|0:48:1:51,51,0|0:49:3:58,50,1|2:21:6:23,27
chr1	100	200	0|0:54:7:56,60,0/1:35:4
chr2	0	100	NULL
chr2	100	200	0|0:54:7:56,60,0/1:35:4
chr3	0	100	NULL
chr3	100	200	NULL" > exp
$BT map -a ivls.bed -b test.vcf -null NULL -c 10 -o collapse > obs
check obs exp

###########################################################
#  Test -s
############################################################
echo -e "    map.t34..\c"
echo \
"chr1	0	10	a1	10	-	.
chr2	10	20	a7	2	+	.
chr3	120	130	a9	4	-	." > exp
$BT map -a ivls2.bed -b values.bed -c 4 -o collapse -s > obs
check obs exp

###########################################################
#  Test -S
############################################################
echo -e "    map.t35..\c"
echo \
"chr1	0	10	a1	10	-	a1
chr2	10	20	a7	2	+	.
chr3	120	130	a9	4	-	a8" > exp
$BT map -a ivls2.bed -b values.bed -c 4 -o collapse -S > obs
check obs exp

###########################################################
#  Test -f 0.1
############################################################
echo -e "    map.t36..\c"
echo \
"chr1	0	10	a1	10	-	a1
chr2	10	20	a7	2	+	.
chr3	120	130	a9	4	-	a8" > exp
$BT map -a ivls2.bed -b values5.bed -c 4 -o collapse -f 0.1 > obs
check obs exp

###########################################################
#  Test -f 0.7
############################################################
echo -e "    map.t37..\c"
echo \
"chr1	0	10	a1	10	-	.
chr2	10	20	a7	2	+	.
chr3	120	130	a9	4	-	a8" > exp
$BT map -a ivls2.bed -b values5.bed -c 4 -o collapse -f 0.7 > obs
check obs exp

###########################################################
#  Test -f 0.9
############################################################
echo -e "    map.t38..\c"
echo \
"chr1	0	10	a1	10	-	.
chr2	10	20	a7	2	+	.
chr3	120	130	a9	4	-	." > exp
$BT map -a ivls2.bed -b values5.bed -c 4 -o collapse -f 0.9 > obs
check obs exp

###########################################################
#  Test -g
############################################################
echo -e "    map.t39..\c"
echo \
"chr1	10	20	chr1
chr2	10	20	chr2
chr2	200	300	.
chr10	10	20	chr10
chr10	20	30	." > exp
$BT map -a a.vsorted.bed -b b.vsorted.bed -c 1 -o collapse > obs
check obs exp

###########################################################
#  Test -g
############################################################
echo -e "    map.t40..\c"
echo \
"chr1	10	20	chr1
chr2	10	20	chr2
chr2	200	300	.
chr10	10	20	chr10
chr10	20	30	." > exp
$BT map -g genome -a a.vsorted.bed -b b.vsorted.bed -c 1 -o collapse > obs
check obs exp

###########################################################
#  Test invalid column
############################################################
echo -e "    map.t41..\c"
echo \
"
*****
***** ERROR: Requested column 41, but database file test.vcf only has fields 1 - 12." > exp
$BT map -a ivls.bed -b test.vcf -c 41 -o collapse 2>&1 > /dev/null | head -3> obs
check obs exp
rm obs exp

###########################################################
#  Test invalid column
############################################################
echo -e "    map.t42..\c"
echo \
"
*****
***** ERROR: Requested column -1, but database file test.vcf only has fields 1 - 12." > exp
$BT map -a ivls.bed -b test.vcf -c -1 -o collapse 2>&1 > /dev/null | head -3> obs
check obs exp
rm obs exp

###########################################################
#  Test invalid column
############################################################
echo -e "    map.t43..\c"
echo \
"
*****
***** ERROR: Requested column 0, but database file test.vcf only has fields 1 - 12." > exp
$BT map -a ivls.bed -b test.vcf -c 0 -o collapse 2>&1 > /dev/null | head -3> obs
check obs exp
rm obs exp


###########################################################
#
#  DEPRECATED
#  Test that Bam database is not allowed
############################################################
echo -e "    map.t44...\c"
#echo -e "\n*****\n***** ERROR: BAM database file not currently supported for column operations." > exp
#$BT map -a ivls.bed -b values.bam 2> obs
#check obs exp
#rm obs exp
echo ok



###########################################################
#  Test that -split option works correctly
############################################################
echo -e "    map.t45...\c"
echo "chr1	0	50	three_blocks_match	15	+	0	0	0	3	10,10,10,	0,20,40,	." > exp
$BT map -o sum -a three_blocks_match.bed -b three_blocks_nomatch.bed -split > obs
check obs exp
rm obs exp






###########################################################
#
#
#  Tests for multiple columns and operations
#
#
############################################################


###########################################################
#  Test that error is given when ops outnumber columns
############################################################
echo -e "    map.t46...\c"
echo \
"chr1	0	100	3	30
chr1	100	200	1	1
chr2	0	100	0	.
chr2	100	200	0	.
chr3	0	100	3	6
chr3	100	200	1	4" > exp
$BT map -a ivls.bed -b values.bed -o count,sum  > obs
check obs exp
rm obs exp


###########################################################
#  Test that error is given when columns outnumber ops,
# if there are two or more ops.
############################################################
echo -e "    map.t47...\c"
echo \
"
*****
***** ERROR: There are 3 columns given, but there are 2 operations."  > exp
$BT map -a ivls.bed -b values.bed -c 5,1,2 -o count,sum 2>&1 > /dev/null | head -3 > obs
check obs exp
rm obs exp


###########################################################
#  Test that numeric ops for non-numeric columns are
# allowed, but give a warning
############################################################
echo -e "    map.t48...\c"
echo \
" ***** WARNING: Non numeric value chr1 in 1.
 ***** WARNING: Non numeric value chr1 in 1.
 ***** WARNING: Non numeric value chr3 in 1.
 ***** WARNING: Non numeric value chr3 in 1." > exp
$BT map -a ivls.bed -b values.bed -c 1 -o sum 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs exp


###########################################################
#  Test that multiple columns are allowed with a 
# single operation
############################################################
echo -e "    map.t49...\c"
echo \
"chr1	0	100	65	9
chr1	100	200	1	7
chr2	0	100	.	.
chr2	100	200	.	.
chr3	0	100	6	7
chr3	100	200	8	23" > exp
$BT map -a ivls.bed -b values4.bed -c 5,7 -o sum > obs
check obs exp
rm obs exp


###########################################################
#  Test that multiple columns are allowed with an
#  equal number of ops that aren't all the same
############################################################
echo -e "    map.t50...\c"
echo \
"chr1	0	100	13.5	65	9
chr1	100	200	120	1	7
chr2	0	100	.	.	.
chr2	100	200	.	.	.
chr3	0	100	10	6	7
chr3	100	200	120	8	23" > exp
$BT map -a ivls.bed -b values4.bed -c 2,5,7 -o mean,sum,sum > obs
check obs exp
rm obs exp


###########################################################
#  Test stdev
############################################################
echo -e "    map.t51...\c"
echo \
"chr1	0	100	3.593976442
chr1	100	200	0
chr2	0	100	.
chr2	100	200	.
chr3	0	100	8.730533902
chr3	100	200	0.5" > exp
$BT map -a ivls.bed -b values4.bed -c 7 -o stdev > obs
check obs exp
rm obs exp

###########################################################
#  Test sstdev
############################################################
echo -e "    map.t52...\c"
echo \
"chr1	0	100	3.937003937
chr1	100	200	.
chr2	0	100	.
chr2	100	200	.
chr3	0	100	10.69267662
chr3	100	200	0.7071067812" > exp
$BT map -a ivls.bed -b values4.bed -c 7 -o sstdev > obs
check obs exp
rm obs exp

###########################################################
#  Test BAM file as DB
############################################################
echo -e "    map.t53...\c"
echo \
"chr1	10000	12000	2.5
chr1	15000	20000	11.44444444" > exp
$BT map -a d.bed -b fullFields.bam -c 5 -o mean > obs
check exp obs
rm exp obs


###########################################################
#  Bug 262: Test stranded map with BedPlus records. 
# -s (lowercase) should give results
############################################################
echo -e "    map.t54...\c"
echo \
"1	3215742	3216021	.	0	-	1
1	3217007	3218115	.	0	-	1" > exp
$BT map -a bug262_a.bed -b bug262_b.bed -s > obs
check exp obs
rm exp obs


###########################################################
#  Bug 262: Test stranded map with BedPlus records. 
# -S (uppercase) should NOT give results
############################################################
echo -e "    map.t55...\c"
echo \
"1	3215742	3216021	.	0	-	.
1	3217007	3218115	.	0	-	." > exp
$BT map -a bug262_a.bed -b bug262_b.bed -S > obs
check exp obs
rm exp obs



###########################################################
#
# NEW TESTS: operations, formatting, options, record types
#
###########################################################

###########################################################
#  map.new.t001: sum adds decimal values; an A with no overlap
#  gets the null value "."
#  chr1:0-100 overlaps scores 1.5, 2.25 and -0.75: sum 3.
#  chr2:0-100 has no B.
###########################################################
echo -e "    map.new.t001...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1.5\nchr1\t20\t30\tb2\t2.25\nchr1\t30\t40\tb3\t-0.75\n' > b.tmp.bed
printf 'chr1\t0\t100\t3\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t003: count_distinct counts unique strings; 0 for no
#  overlap
#  Names b, a, b have 2 distinct values. pins the current output:
#  count_distinct prints 0, not ".", for no overlap; map.rst says an A with no
#  overlap gets "." or the -null value (bugs.md #16).
###########################################################
echo -e "    map.new.t003...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t1\nchr1\t20\t30\ta\t2\nchr1\t30\t40\tb\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\t2\nchr2\t0\t100\t0\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o count_distinct > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t004: min with negative and decimal values
#  Scores 3, -2.5, 10: min -2.5.
###########################################################
echo -e "    map.new.t004...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t3\nchr1\t20\t30\tb2\t-2.5\nchr1\t30\t40\tb3\t10\n' > b.tmp.bed
printf 'chr1\t0\t100\t-2.5\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o min > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t005: max compares numerically
#  Scores 3, -2.5, 10: max 10 (a string comparison would pick 3).
###########################################################
echo -e "    map.new.t005...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t3\nchr1\t20\t30\tb2\t-2.5\nchr1\t30\t40\tb3\t10\n' > b.tmp.bed
printf 'chr1\t0\t100\t10\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o max > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t006: absmin keeps the fractional part
#  Scores 3, -2.5, 10: absolute values 3, 2.5, 10; min 2.5.
###########################################################
echo -e "    map.new.t006...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t3\nchr1\t20\t30\tb2\t-2.5\nchr1\t30\t40\tb3\t10\n' > b.tmp.bed
printf 'chr1\t0\t100\t2.5\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o absmin > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t007: absmax reports the absolute value
#  Scores -12, 3, 10: absolute values 12, 3, 10; max 12.
###########################################################
echo -e "    map.new.t007...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t-12\nchr1\t20\t30\tb2\t3\nchr1\t30\t40\tb3\t10\n' > b.tmp.bed
printf 'chr1\t0\t100\t12\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o absmax > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t008: mean prints 10 significant digits by default
#  Scores 1, 2, 4: mean 7/3 = 2.3333333333... This is the one
#  test of the default precision. pins the current output: the
#  default precision is 10 significant digits; the usage text says
#  5 (bugs.md #10).
###########################################################
echo -e "    map.new.t008...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\nchr1\t30\t40\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t2.333333333\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o mean > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t009: median of an odd count sorts numerically
#  Scores 10, 9, 2 sort to 2, 9, 10: median 9 (a string sort,
#  "10" "2" "9", would give 2).
###########################################################
echo -e "    map.new.t009...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t10\nchr1\t20\t30\tb2\t9\nchr1\t30\t40\tb3\t2\n' > b.tmp.bed
printf 'chr1\t0\t100\t9\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o median > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t010: median of an even count averages the middle two
#  Scores 4, 1, 3, 10 sort to 1, 3, 4, 10: (3 + 4) / 2 = 3.5.
###########################################################
echo -e "    map.new.t010...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t4\nchr1\t20\t30\tb2\t1\nchr1\t30\t40\tb3\t3\nchr1\t40\t50\tb4\t10\n' > b.tmp.bed
printf 'chr1\t0\t100\t3.5\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o median > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t011: mode breaks a frequency tie with the smallest
#  string
#  Names b, a, b, a, c: a and b both occur twice. pins the current
#  output: the lexicographically smallest of the tied values, a.
###########################################################
echo -e "    map.new.t011...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t1\nchr1\t20\t30\ta\t1\nchr1\t30\t40\tb\t1\nchr1\t40\t50\ta\t1\nchr1\t50\t60\tc\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\ta\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o mode > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t012: mode compares numbers as strings
#  Scores 9, 10, 9, 10 tie; "10" sorts before "9" as a string.
#  pins the current output: 10.
###########################################################
echo -e "    map.new.t012...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t9\nchr1\t20\t30\tb2\t10\nchr1\t30\t40\tb3\t9\nchr1\t40\t50\tb4\t10\n' > b.tmp.bed
printf 'chr1\t0\t100\t10\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o mode > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t013: antimode picks the least frequent value, smallest
#  string first on ties
#  Names b, c, a, b: a and c occur once, b twice. pins the current
#  output: a.
###########################################################
echo -e "    map.new.t013...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t1\nchr1\t20\t30\tc\t1\nchr1\t30\t40\ta\t1\nchr1\t40\t50\tb\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\ta\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o antimode > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t014: collapse lists values in B order with duplicates
#  Names b, a, b -> b,a,b.
###########################################################
echo -e "    map.new.t014...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t1\nchr1\t20\t30\ta\t2\nchr1\t30\t40\tb\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\tb,a,b\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o collapse > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t015: distinct lists unique values sorted as strings
#  Scores 10, 9, 2, 9 -> unique 10, 2, 9 in string order. pins the
#  current output: the usage text does not give the order.
###########################################################
echo -e "    map.new.t015...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t10\nchr1\t20\t30\tb2\t9\nchr1\t30\t40\tb3\t2\nchr1\t40\t50\tb4\t9\n' > b.tmp.bed
printf 'chr1\t0\t100\t10,2,9\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o distinct > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t016: distinct_only keeps only values seen once
#  Names a, b, b, c: a and c occur once.
###########################################################
echo -e "    map.new.t016...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\ta\t1\nchr1\t20\t30\tb\t2\nchr1\t30\t40\tb\t3\nchr1\t40\t50\tc\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\ta,c\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o distinct_only > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t017: first reports the value of the first B in file
#  order
#  Names x, y, z -> x.
###########################################################
echo -e "    map.new.t017...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tx\t1\nchr1\t20\t30\ty\t2\nchr1\t30\t40\tz\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\tx\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o first > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t018: last reports the value of the last B in file order
#  Names x, y, z -> z.
###########################################################
echo -e "    map.new.t018...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tx\t1\nchr1\t20\t30\ty\t2\nchr1\t30\t40\tz\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\tz\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o last > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t019: stdev is the population standard deviation
#  Scores 2,4,4,4,5,5,7,9: mean 5, squared deviations sum to 32,
#  sqrt(32/8) = 2.
###########################################################
echo -e "    map.new.t019...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t1\t2\tb\t2\nchr1\t2\t3\tb\t4\nchr1\t3\t4\tb\t4\nchr1\t4\t5\tb\t4\nchr1\t5\t6\tb\t5\nchr1\t6\t7\tb\t5\nchr1\t7\t8\tb\t7\nchr1\t8\t9\tb\t9\n' > b.tmp.bed
printf 'chr1\t0\t100\t2\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o stdev > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t020: sstdev is the sample standard deviation
#  Same scores as t019: sqrt(32/7) = 2.13808993529..., printed with
#  -prec 10 significant digits.
###########################################################
echo -e "    map.new.t020...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t1\t2\tb\t2\nchr1\t2\t3\tb\t4\nchr1\t3\t4\tb\t4\nchr1\t4\t5\tb\t4\nchr1\t5\t6\tb\t5\nchr1\t6\t7\tb\t5\nchr1\t7\t8\tb\t7\nchr1\t8\t9\tb\t9\n' > b.tmp.bed
printf 'chr1\t0\t100\t2.138089935\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sstdev -prec 10 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t021: stdev and sstdev of a single value
#  Population: sqrt(0/1) = 0. Sample: 0/0 is undefined, so the
#  null value is printed.
###########################################################
echo -e "    map.new.t021...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t6\n' > b.tmp.bed
printf 'chr1\t0\t100\t0\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o stdev,sstdev > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t022: freqasc lists value:count by increasing count
#  Names b, a, b, a, c: c:1, then a:2 and b:2. pins the current
#  output: values with equal counts are in string order.
###########################################################
echo -e "    map.new.t022...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t1\nchr1\t20\t30\ta\t1\nchr1\t30\t40\tb\t1\nchr1\t40\t50\ta\t1\nchr1\t50\t60\tc\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\tc:1,a:2,b:2\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o freqasc > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t023: freqdesc lists value:count by decreasing count
#  Same names as t022: a:2, b:2, then c:1. pins the current output:
#  values with equal counts are in string order.
###########################################################
echo -e "    map.new.t023...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t1\nchr1\t20\t30\ta\t1\nchr1\t30\t40\tb\t1\nchr1\t40\t50\ta\t1\nchr1\t50\t60\tc\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\ta:2,b:2,c:1\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o freqdesc > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t024: concat joins values in B order with no delimiter
#  Names a, b, c -> abc. concat is accepted by the parser though
#  the usage text does not list it. pins the current output: concat
#  uses no delimiter; map.rst says "print a comma separated list"
#  (bugs.md #15).
###########################################################
echo -e "    map.new.t024...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\ta\t1\nchr1\t20\t30\tb\t2\nchr1\t30\t40\tc\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\tabc\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o concat > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t025: distinct_sort_num sorts unique values numerically
#  Scores 10, 9, 2, 9 -> 2,9,10.
###########################################################
echo -e "    map.new.t025...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t10\nchr1\t20\t30\tb2\t9\nchr1\t30\t40\tb3\t2\nchr1\t40\t50\tb4\t9\n' > b.tmp.bed
printf 'chr1\t0\t100\t2,9,10\nchr2\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o distinct_sort_num > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t026: distinct_sort_num_desc sorts unique values
#  numerically, descending
#  Scores 10, 9, 2, 9 -> 10,9,2.
###########################################################
echo -e "    map.new.t026...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t10\nchr1\t20\t30\tb2\t9\nchr1\t30\t40\tb3\t2\nchr1\t40\t50\tb4\t9\n' > b.tmp.bed
printf 'chr1\t0\t100\t10,9,2\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o distinct_sort_num_desc > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t027: -delim sets the collapse delimiter
#  Names b, a, b with -delim "|" -> b|a|b.
###########################################################
echo -e "    map.new.t027...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t1\nchr1\t20\t30\ta\t2\nchr1\t30\t40\tb\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\tb|a|b\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o collapse -delim "|" > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t028: -delim also separates distinct, freqdesc and
#  distinct_sort_num items
#  Names b, a, b and scores 3, 1, 3: distinct a;b, freqdesc
#  b:2;a:1, distinct_sort_num 1;3.
###########################################################
echo -e "    map.new.t028...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t3\nchr1\t20\t30\ta\t1\nchr1\t30\t40\tb\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\ta;b\tb:2;a:1\t1;3\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4,4,5 -o distinct,freqdesc,distinct_sort_num -delim ";" > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t029: -null replaces "." for numeric and list operations
#  but count and count_distinct still print 0
#  chr1:0-100 overlaps b (1) and a (2); chr2:0-100 overlaps nothing.
#  pins the current output: count and count_distinct print 0 and
#  ignore -null; map.rst says an A with no
#  overlap gets "." or the -null value (bugs.md #16).
###########################################################
echo -e "    map.new.t029...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb\t1\nchr1\t20\t30\ta\t2\n' > b.tmp.bed
printf 'chr1\t0\t100\t3\tb,a\t2\t2\nchr2\t0\t100\tNA\tNA\t0\t0\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5,4,5,4 -o sum,collapse,count,count_distinct -null NA > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t030: -prec sets significant digits
#  Scores 1, 2, 4: mean 7/3 with -prec 3 -> 2.33. pins the current
#  output: -prec counts significant digits; the usage text says
#  "decimal precision" (bugs.md #10).
###########################################################
echo -e "    map.new.t030...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\nchr1\t30\t40\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t2.33\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o mean -prec 3 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t031: a sum with more digits than -prec switches to
#  scientific notation
#  Scores 12000 and 345: sum 12345 with -prec 2. pins the current
#  output: 1.2e+04.
###########################################################
echo -e "    map.new.t031...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t12000\nchr1\t20\t30\tb2\t345\n' > b.tmp.bed
printf 'chr1\t0\t100\t1.2e+04\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum -prec 2 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t032: an 11-digit sum exceeds -prec 10
#  Scores 12345678900 and 1: sum 12345678901 has 11 digits. pins
#  the current output: 1.23456789e+10.
###########################################################
echo -e "    map.new.t032...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t12345678900\nchr1\t20\t30\tb2\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\t1.23456789e+10\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum -prec 10 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t033: a 10-digit sum prints every digit at -prec 10
#  Scores 1234567800 and 90: sum 1234567890.
###########################################################
echo -e "    map.new.t033...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1234567800\nchr1\t20\t30\tb2\t90\n' > b.tmp.bed
printf 'chr1\t0\t100\t1234567890\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum -prec 10 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t034: one column with several operations, output in the
#  order the operations are given
#  Scores 1, 2, 4 with -prec 10: sum 7, count 3, mean 2.333333333,
#  min 1, max 4.
###########################################################
echo -e "    map.new.t034...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\nchr1\t30\t40\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t7\t3\t2.333333333\t1\t4\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum,count,mean,min,max -prec 10 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t035: several columns with one operation
#  -c 2,3 -o max: largest B start 30 and largest B end 40.
###########################################################
echo -e "    map.new.t035...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\nchr1\t30\t40\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t30\t40\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 2,3 -o max > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t036: column/operation pairs are applied in order
#  -c 5,4 -o sum,collapse: 1+2+4 = 7, then names b1,b2,b3.
###########################################################
echo -e "    map.new.t036...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\nchr1\t30\t40\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t7\tb1,b2,b3\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5,4 -o sum,collapse > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t037: -ops is an alias for -o
#  Same input as t036 with -c 5 -ops max: 4.
###########################################################
echo -e "    map.new.t037...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\nchr1\t30\t40\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t4\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -ops max > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t038: -prec 0 is an error
###########################################################
echo -e "    map.new.t038...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\n' > b.tmp.bed
printf '***** ERROR: -prec must be followed by a positive integer. Exiting. *****\nexit 1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -o mean -prec 0 > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t039: more than one -b file is an error
###########################################################
echo -e "    map.new.t039...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\n' > b1.tmp.bed
printf 'chr1\t10\t20\tb2\t2\n' > b2.tmp.bed
printf '***** ERROR: multiple database files currently not supported for map. *****\nexit 1\n' > exp
$BT map -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  map.new.t040: -e with -F but no -f is an error
###########################################################
echo -e "    map.new.t040...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\n' > b.tmp.bed
printf '***** ERROR: -e must be used with -f. *****\nexit 1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -F 0.5 -e > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t041: -r with -F is an error
###########################################################
echo -e "    map.new.t041...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\n' > b.tmp.bed
printf '***** ERROR: -r must be used solely with -f. *****\nexit 1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -F 0.5 -r > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t042: -s uses only B on the same known strand
#  A is +. b1 (+, 1), b2 (-, 2), b3 (+, 4), b4 (., 8): -s sums
#  1 + 4 = 5.
###########################################################
echo -e "    map.new.t042...\c"
printf 'chr1\t0\t100\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\t+\nchr1\t20\t30\tb2\t2\t-\nchr1\t30\t40\tb3\t4\t+\nchr1\t40\t50\tb4\t8\t.\n' > b.tmp.bed
printf 'chr1\t0\t100\ta1\t0\t+\t5\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -s > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t043: -S uses only B on the opposite known strand
#  Same input as t042: only b2 (-) qualifies: 2. The unknown-strand
#  b4 is excluded.
###########################################################
echo -e "    map.new.t043...\c"
printf 'chr1\t0\t100\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\t+\nchr1\t20\t30\tb2\t2\t-\nchr1\t30\t40\tb3\t4\t+\nchr1\t40\t50\tb4\t8\t.\n' > b.tmp.bed
printf 'chr1\t0\t100\ta1\t0\t+\t2\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -S > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t044: -f keeps B covering exactly the fraction of A
#  A chr1:0-100. b1 chr1:50-150 overlaps 50 bp = 0.50 (kept at
#  -f 0.5), b2 chr1:60-200 overlaps 40 bp = 0.40 (dropped): sum 1.
###########################################################
echo -e "    map.new.t044...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\tb1\t1\nchr1\t60\t200\tb2\t2\n' > b.tmp.bed
printf 'chr1\t0\t100\t1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -f 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t045: -F tests the fraction of each B
#  A chr1:0-100. b1 chr1:50-150: 50/100 = 0.50 of B (kept);
#  b2 chr1:80-200: 20/120 = 0.17 (dropped); b3 chr1:90-100:
#  10/10 = 1.0 (kept). -F 0.5 sums 1 + 4 = 5.
###########################################################
echo -e "    map.new.t045...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\tb1\t1\nchr1\t80\t200\tb2\t2\nchr1\t90\t100\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t5\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -F 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t046: -f with -r requires the fraction of A and of B
#  Same input as t045 with -f 0.5 -r: b1 0.50/0.50 (kept), b2
#  0.20/0.17, b3 0.10/1.0 (dropped): sum 1.
###########################################################
echo -e "    map.new.t046...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\tb1\t1\nchr1\t80\t200\tb2\t2\nchr1\t90\t100\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -f 0.5 -r > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t047: -e accepts B passing either fraction
#  Same input as t045 with -f 0.6 -F 0.6 -e: b1 0.50/0.50 fails
#  both, b2 0.20/0.17 fails both, b3 0.10/1.0 passes -F: sum 4.
###########################################################
echo -e "    map.new.t047...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t50\t150\tb1\t1\nchr1\t80\t200\tb2\t2\nchr1\t90\t100\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t4\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -f 0.6 -F 0.6 -e > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t048: a BED12 A without -split uses its whole span
#  A chr1:0-100 has blocks 0-10 and 90-100. b1 chr1:40-50 (in the
#  gap, score 1) and b2 chr1:5-15 (in block 1, score 2) both
#  overlap the span: sum 3. All 12 fields of A are printed.
###########################################################
echo -e "    map.new.t048...\c"
printf 'chr1\t0\t100\ta1\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,\n' > a.tmp.bed
printf 'chr1\t5\t15\tb2\t2\nchr1\t40\t50\tb1\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\ta1\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,\t3\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t049: -split with a BED12 A drops B in the gap
#  Same input as t048: only b2 overlaps a block: sum 2.
###########################################################
echo -e "    map.new.t049...\c"
printf 'chr1\t0\t100\ta1\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,\n' > a.tmp.bed
printf 'chr1\t5\t15\tb2\t2\nchr1\t40\t50\tb1\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\ta1\t0\t+\t0\t100\t0\t2\t10,10,\t0,90,\t2\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -split > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t050: -split with a BED12 B drops it when A lies in its
#  gap
#  B chr1:0-100 has blocks 0-10 and 90-100 (score 7). A1
#  chr1:40-50 lies in the gap: null. A2 chr1:95-120 overlaps the
#  second block: 7.
###########################################################
echo -e "    map.new.t050...\c"
printf 'chr1\t40\t50\nchr1\t95\t120\n' > a.tmp.bed
printf 'chr1\t0\t100\tb1\t7\t+\t0\t100\t0\t2\t10,10,\t0,90,\n' > b.tmp.bed
printf 'chr1\t40\t50\t.\nchr1\t95\t120\t7\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -split > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t051: -header prints the query's "#" and "track" lines
###########################################################
echo -e "    map.new.t051...\c"
printf '# comment\ntrack name=a\nchr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\n' > b.tmp.bed
printf '# comment\ntrack name=a\nchr1\t0\t100\t1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t052: a non-numeric value in a numeric operation gives
#  the null value and a warning
#  Scores 1, x, 3: sum is undefined, so "." on stdout; stderr has
#  the warning naming x and column 5. The test appends stderr
#  after stdout.
###########################################################
echo -e "    map.new.t052...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\tx\nchr1\t30\t40\tb3\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\t.\n ***** WARNING: Non numeric value x in 5.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum > obs 2> err.tmp
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t053: one warning per A record, naming the last
#  non-numeric value
#  Scores x, 2, y for one A: one warning, for y. pins the current
#  output.
###########################################################
echo -e "    map.new.t053...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\tx\nchr1\t20\t30\tb2\t2\nchr1\t30\t40\tb3\ty\n' > b.tmp.bed
printf 'chr1\t0\t100\t.\n ***** WARNING: Non numeric value y in 5.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum > obs 2> err.tmp
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t054: non-numeric operations on a text column give no
#  warning
#  count and collapse of names b1, b2: 2 and b1,b2; stderr empty.
###########################################################
echo -e "    map.new.t054...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\n' > b.tmp.bed
printf 'chr1\t0\t100\t2\tb1,b2\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o count,collapse > obs 2> err.tmp
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t055: exponent notation is numeric
#  Scores 1e2 and 2.5E1: sum 100 + 25 = 125.
###########################################################
echo -e "    map.new.t055...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1e2\nchr1\t20\t30\tb2\t2.5E1\n' > b.tmp.bed
printf 'chr1\t0\t100\t125\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t056: a GFF query prints as read, with 1-based start
#  converted for the overlap test
#  GFF chr1 1-100 is 0-based chr1:0-100. b0 chr1:0-1 (1) and b1
#  chr1:99-110 (5) overlap it; b2 chr1:100-110 (7) is book-ended.
#  Sum 6 (reading 1 as 0-based, A would be book-ended to b0: 5).
###########################################################
echo -e "    map.new.t056...\c"
printf 'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=a1\n' > a.tmp.gff
printf 'chr1\t0\t1\tb0\t1\nchr1\t99\t110\tb1\t5\nchr1\t100\t110\tb2\t7\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=a1\t6\n' > exp
$BT map -a a.tmp.gff -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  map.new.t057: a VCF query prints as read, without its header
#  POS 101 REF A is 0-based chr1:100-101. b2 chr1:100-110 (3)
#  overlaps it; b1 chr1:101-110 (4) starts at its end and does not.
#  Sum 3 (reading POS as 0-based, A chr1:101-102 would overlap both:
#  7).
###########################################################
echo -e "    map.new.t057...\c"
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\trs1\tA\tG\t50\tPASS\tDP=3\n' > a.tmp.vcf
printf 'chr1\t100\t110\tb2\t3\nchr1\t101\t110\tb1\t4\n' > b.tmp.bed
printf 'chr1\t101\trs1\tA\tG\t50\tPASS\tDP=3\t3\n' > exp
$BT map -a a.tmp.vcf -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  map.new.t058: a BED6+3 query prints every column before the
#  value
###########################################################
echo -e "    map.new.t058...\c"
printf 'chr1\t0\t100\ta1\t9\t-\tp\tq\tr\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\n' > b.tmp.bed
printf 'chr1\t0\t100\ta1\t9\t-\tp\tq\tr\t3\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t059: a zero-length query inside B gets its value and is
#  printed with its original coordinates
#  A chr1:50-50 lies inside B chr1:10-100 (score 6).
###########################################################
echo -e "    map.new.t059...\c"
printf 'chr1\t50\t50\n' > a.tmp.bed
printf 'chr1\t10\t100\tb1\t6\n' > b.tmp.bed
printf 'chr1\t50\t50\t6\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t060: book-ended B records are not overlaps
#  A chr1:100-200. chr1:50-100 (1) and chr1:200-250 (2) only touch
#  A; chr1:150-160 (4) overlaps: sum 4.
###########################################################
echo -e "    map.new.t060...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t100\tb1\t1\nchr1\t150\t160\tb3\t4\nchr1\t200\t250\tb2\t2\n' > b.tmp.bed
printf 'chr1\t100\t200\t4\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t061: a B spanning two A records counts for both
#  B chr1:50-150 (3) overlaps chr1:0-100 and chr1:100-200.
###########################################################
echo -e "    map.new.t061...\c"
printf 'chr1\t0\t100\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t150\tb1\t3\n' > b.tmp.bed
printf 'chr1\t0\t100\t3\nchr1\t100\t200\t3\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t062: identical A records each get the value
#  Two copies of chr1:0-100 overlap B scores 1 and 2.
###########################################################
echo -e "    map.new.t062...\c"
printf 'chr1\t0\t100\nchr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr1\t20\t30\tb2\t2\n' > b.tmp.bed
printf 'chr1\t0\t100\t3\nchr1\t0\t100\t3\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t063: B chromosomes absent from A are skipped
#  A has chr1 and chr3; B also has chr2 (score 2) between them.
#  chr1 gets 1, chr3 gets 4.
###########################################################
echo -e "    map.new.t063...\c"
printf 'chr1\t0\t100\nchr3\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\nchr2\t10\t20\tb2\t2\nchr3\t10\t20\tb3\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t1\nchr3\t0\t100\t4\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t064: GFF database, score column, 1-based start
#  GFF e1 101-150 is 0-based chr1:100-150, e2 151-250 is
#  chr1:150-250. A chr1:0-101 overlaps e1 by 1 bp: 2.5 (reading 101
#  as 0-based, e1 would be book-ended: "."). A chr1:100-200 overlaps
#  both: 2.5 + 1.5 = 4.
###########################################################
echo -e "    map.new.t064...\c"
printf 'chr1\t0\t101\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\tsrc\texon\t101\t150\t2.5\t+\t.\tID=e1\nchr1\tsrc\texon\t151\t250\t1.5\t+\t.\tID=e2\n' > b.tmp.gff
printf 'chr1\t0\t101\t2.5\nchr1\t100\t200\t4\n' > exp
$BT map -a a.tmp.bed -b b.tmp.gff -c 6 -o sum > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  map.new.t065: VCF database, QUAL column, 1-based POS
#  POS 100 is 0-based chr1:99-100 (QUAL 10), POS 101 is
#  chr1:100-101 (QUAL 20). A chr1:0-100 gets 10, chr1:100-200
#  gets 20.
###########################################################
echo -e "    map.new.t065...\c"
printf 'chr1\t0\t100\nchr1\t100\t200\n' > a.tmp.bed
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\trs1\tA\tG\t10\tPASS\t.\nchr1\t101\trs2\tC\tT\t20\tPASS\t.\n' > b.tmp.vcf
printf 'chr1\t0\t100\t10\nchr1\t100\t200\t20\n' > exp
$BT map -a a.tmp.bed -b b.tmp.vcf -c 6 -o sum > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  map.new.t066: bedGraph database, value column
#  Values 1 and 4 over chr1:0-100: mean 2.5.
###########################################################
echo -e "    map.new.t066...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t0\t50\t1\nchr1\t50\t100\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t2.5\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 4 -o mean > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t067: with -g, inputs follow the genome file's
#  chromosome order
#  The genome lists chr2 before chr1, and both files follow it.
###########################################################
echo -e "    map.new.t067...\c"
printf 'chr2\t1000\nchr1\t1000\n' > g.tmp
printf 'chr2\t0\t100\nchr1\t0\t100\n' > a.tmp.bed
printf 'chr2\t10\t20\tb1\t1\nchr1\t10\t20\tb2\t2\n' > b.tmp.bed
printf 'chr2\t0\t100\t1\nchr1\t0\t100\t2\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t068: with -g, a file in another chromosome order is an
#  error
#  The genome lists chr2 before chr1, but A has chr1 first.
###########################################################
echo -e "    map.new.t068...\c"
printf 'chr2\t1000\nchr1\t1000\n' > g.tmp
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr2\t10\t20\tb1\t1\nchr1\t10\t20\tb2\t2\n' > b.tmp.bed
printf 'Error: Sorted input specified, but the file a.tmp.bed has the following record with a different sort order than the genomeFile g.tmp\nchr2\t0\t100\nexit 1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -g g.tmp > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp g.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t069: an unsorted database is an error
###########################################################
echo -e "    map.new.t069...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t30\t40\tb1\t1\nchr1\t10\t20\tb2\t2\n' > b.tmp.bed
printf 'Error: Sorted input specified, but the file b.tmp.bed has the following out of order record\nchr1\t10\t20\tb2\t2\nexit 1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t070: the query can be read from stdin
#  Same records as t001.
###########################################################
echo -e "    map.new.t070...\c"
printf 'chr1\t0\t100\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1.5\nchr1\t20\t30\tb2\t2.25\nchr1\t30\t40\tb3\t-0.75\n' > b.tmp.bed
printf 'chr1\t0\t100\t3\nchr2\t0\t100\t.\n' > exp
cat a.tmp.bed | $BT map -a - -b b.tmp.bed -c 5 -o sum > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t071: a BAM query with -bed is printed as BED12
#  values.bam holds a1-a4 on chr1 at POS 1, 11, 21, 121 and a5-a8
#  on chr3 at the same positions, each 10M, MAPQ 255, flag 0: so
#  0-based chr1:0-10, 10-20, 20-30, 120-130 (likewise chr3), name,
#  score 255, strand +. B chr1:5-15 (2) overlaps a1 and a2; B
#  chr3:0-5 (7) overlaps a5. pins the current output: the key is
#  printed as BED12 (thick start/end = start/end, itemRgb 0,0,0,
#  one block of the aligned length).
###########################################################
echo -e "    map.new.t071...\c"
printf 'chr1\t5\t15\tb1\t2\nchr3\t0\t5\tb2\t7\n' > b.tmp.bed
printf '%b\n' \
	'chr1\t0\t10\ta1\t255\t+\t0\t10\t0,0,0\t1\t10,\t0,\t2' \
	'chr1\t10\t20\ta2\t255\t+\t10\t20\t0,0,0\t1\t10,\t0,\t2' \
	'chr1\t20\t30\ta3\t255\t+\t20\t30\t0,0,0\t1\t10,\t0,\t.' \
	'chr1\t120\t130\ta4\t255\t+\t120\t130\t0,0,0\t1\t10,\t0,\t.' \
	'chr3\t0\t10\ta5\t255\t+\t0\t10\t0,0,0\t1\t10,\t0,\t7' \
	'chr3\t10\t20\ta6\t255\t+\t10\t20\t0,0,0\t1\t10,\t0,\t.' \
	'chr3\t20\t30\ta7\t255\t+\t20\t30\t0,0,0\t1\t10,\t0,\t.' \
	'chr3\t120\t130\ta8\t255\t+\t120\t130\t0,0,0\t1\t10,\t0,\t.' > exp
$BT map -a values.bam -b b.tmp.bed -bed -c 5 -o sum > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  map.new.t072: an input buffer shorter than one line gives the
#  same result
#  -iobuf 8 forces many refills per line. Same records as t001
#  plus a second chr1 A and a chr2 B: chr1:0-100 gets 1.5 + 2.25 =
#  3.75 (collapse 1.5,2.25), chr1:100-200 gets -0.75, chr2:0-100
#  gets 4.
###########################################################
echo -e "    map.new.t072...\c"
printf 'chr1\t0\t100\nchr1\t100\t200\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1.5\nchr1\t20\t30\tb2\t2.25\nchr1\t150\t160\tb3\t-0.75\nchr2\t50\t60\tb4\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t3.75\t1.5,2.25\nchr1\t100\t200\t-0.75\t-0.75\nchr2\t0\t100\t4\t4\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum,collapse -iobuf 8 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t073: -nobuf writes the same lines
#  Same input and result as t072, without -iobuf.
###########################################################
echo -e "    map.new.t073...\c"
printf 'chr1\t0\t100\nchr1\t100\t200\nchr2\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1.5\nchr1\t20\t30\tb2\t2.25\nchr1\t150\t160\tb3\t-0.75\nchr2\t50\t60\tb4\t4\n' > b.tmp.bed
printf 'chr1\t0\t100\t3.75\t1.5,2.25\nchr1\t100\t200\t-0.75\t-0.75\nchr2\t0\t100\t4\t4\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -c 5 -o sum,collapse -nobuf > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t074: an -iobuf below the 8-byte minimum is an error
###########################################################
echo -e "    map.new.t074...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf 'chr1\t10\t20\tb1\t1\n' > b.tmp.bed
printf '***** ERROR: specified buffer size is too small. *****\nexit 1\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -iobuf 4 > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t075: chromosome names with and without "chr" give a
#  warning
#  A uses chr1, B uses 1, so nothing overlaps (null on stdout) and
#  stderr warns about the B record. pins the current output: the
#  warning is printed when found and again at exit. The test
#  appends stderr after stdout.
###########################################################
echo -e "    map.new.t075...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf '1\t10\t20\tb1\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\t.\n***** WARNING: File b.tmp.bed has inconsistent naming convention for record:\n1\t10\t20\tb1\t1\n\n***** WARNING: File b.tmp.bed has inconsistent naming convention for record:\n1\t10\t20\tb1\t1\n\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed > obs 2> err.tmp
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  map.new.t076: -nonamecheck suppresses the naming warning
#  Same input as t075: only the null line, nothing on stderr.
###########################################################
echo -e "    map.new.t076...\c"
printf 'chr1\t0\t100\n' > a.tmp.bed
printf '1\t10\t20\tb1\t1\n' > b.tmp.bed
printf 'chr1\t0\t100\t.\n' > exp
$BT map -a a.tmp.bed -b b.tmp.bed -nonamecheck > obs 2> err.tmp
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

[[ $FAILURES -eq 0 ]] || exit 1;
