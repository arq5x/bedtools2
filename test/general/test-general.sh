BT=${BT-../../bin/bedtools}

FAILURES=0;

check()
{
	if diff $1 $2; then
    	echo ok
		return 1
	else
    	FAILURES=$(expr $FAILURES + 1);
		echo fail
		return 0
	fi
}

###########################################################
#  Enforce non-negative coordinates
###########################################################
echo -e "    general.t01...\c"
echo \
"chr1	1	10
chr1	-1	10" | $BT merge -i - 2> obs
echo "Error: Invalid record in file -. Record is 
chr1	-1	10" > exp
check obs exp
rm obs exp

###########################################################
#  Enforce start <= end
###########################################################
echo -e "    general.t02...\c"
echo \
"chr1	1	2
chr1	10	5" | $BT merge -i - 2> obs
echo "Error: Invalid record in file -. Record is 
chr1	10	5" > exp
check obs exp
rm obs exp

###########################################################
#  Enforce integer coordinates
###########################################################
echo -e "    general.t03...\c"
echo \
"ERROR: file - has non positional records, which are only valid for the groupBy tool." > exp
echo "chr1	.	2" | $BT merge -i - 2> o
head -n 1 o > obs
check obs exp
rm obs exp



###########################################################
#  Enforce tab-separated files
###########################################################
echo -e "    general.t05...\c"
echo \
"chr1 1 2" | $BT merge -i - 2> obs
echo \
"ERROR: file - has non positional records, which are only valid for the groupBy tool." > exp
head -n 1 o > obs
check obs exp
rm obs exp


###########################################################
#  Fail on non-existent files.
###########################################################
echo -e "    general.t06...\c"
$BT merge -i idontexist.bed 2> obs
echo "Error: Unable to open file idontexist.bed. Exiting." > exp
check obs exp
rm obs exp


###########################################################
#  Don't fail on existent, yet empty files.
###########################################################
echo -e "    general.t07...\c"
$BT merge -i empty.bed 2> obs
touch exp
check obs exp
rm obs exp


###########################################################
#  Process gzipped files.
###########################################################
echo -e "    general.t08...\c"
$BT merge -i non-empty.bed.gz > obs
echo "chr1	10	21" > exp
check obs exp
rm obs exp


###########################################################
#  Test GZIP, non-GZIP, and empty file functionality.
###########################################################
echo "chr1	1	100" > a.bed
echo "chr1	75	200" > b.bed
echo -n "" > c.bed
echo "chr1	1	5000" > genome.txt
gzip -c a.bed > a.bed.gz
gzip -c b.bed > b.bed.gz
gzip -c c.bed > c.bed.gz


echo -e "    general.t09...\c"
$BT intersect -a a.bed -b b.bed > obs
echo "chr1	75	100" > exp
check obs exp

echo -e "    general.t10...\c"
$BT intersect -a a.bed.gz -b b.bed.gz > obs
echo "chr1	75	100" > exp
check obs exp

echo -e "    general.t11...\c"
$BT intersect -a a.bed -b b.bed.gz > obs
echo "chr1	75	100" > exp
check obs exp

echo -e "    general.t12...\c"
$BT intersect -a a.bed.gz -b b.bed > obs
echo "chr1	75	100" > exp
check obs exp

echo -e "    general.t13...\c"
$BT intersect -a c.bed -b b.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.t14...\c"
$BT intersect -a c.bed.gz -b b.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.15...\c"
$BT intersect -a c.bed.gz -b c.bed.gz > obs
echo -n "" > exp
check obs exp


echo -e "    general.t16...\c"
$BT subtract -a a.bed -b b.bed > obs
echo "chr1	1	75" > exp
check obs exp

echo -e "    general.t17...\c"
$BT subtract -a a.bed.gz -b b.bed.gz > obs
echo "chr1	1	75" > exp
check obs exp

echo -e "    general.t18...\c"
$BT subtract -a a.bed -b b.bed.gz > obs
echo "chr1	1	75" > exp
check obs exp

echo -e "    general.t19...\c"
$BT subtract -a a.bed.gz -b b.bed > obs
echo "chr1	1	75" > exp
check obs exp

echo -e "    general.t20...\c"
$BT subtract -a c.bed -b b.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.t21...\c"
$BT subtract -a c.bed.gz -b b.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.22...\c"
$BT subtract -a c.bed.gz -b c.bed.gz > obs
echo -n "" > exp
check obs exp


echo -e "    general.t23...\c"
$BT window -a a.bed -b b.bed > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo -e "    general.t24...\c"
$BT window -a a.bed.gz -b b.bed.gz > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo -e "    general.t25...\c"
$BT window -a a.bed -b b.bed.gz > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo -e "    general.t26...\c"
$BT window -a a.bed.gz -b b.bed > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo -e "    general.t27...\c"
$BT window -a c.bed -b b.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.t28...\c"
$BT window -a c.bed.gz -b b.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.29...\c"
$BT window -a c.bed.gz -b c.bed.gz > obs
echo -n "" > exp
check obs exp


echo -e "    general.t30...\c"
$BT closest -a a.bed -b b.bed > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo -e "    general.t31...\c"
$BT closest -a a.bed.gz -b b.bed.gz > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo -e "    general.t32...\c"
$BT closest -a a.bed -b b.bed.gz > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo -e "    general.t33...\c"
$BT closest -a a.bed.gz -b b.bed > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo -e "    general.t34...\c"
$BT closest -a c.bed -b b.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.t35...\c"
$BT closest -a c.bed.gz -b b.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.36...\c"
$BT closest -a c.bed.gz -b c.bed.gz > obs
echo -n "" > exp
check obs exp


echo -e "    general.t37...\c"
$BT merge -i a.bed > obs
echo "chr1	1	100" > exp
check obs exp

echo -e "    general.t38...\c"
$BT merge -i a.bed.gz > obs
echo "chr1	1	100" > exp
check obs exp

echo -e "    general.t39...\c"
$BT merge -i b.bed > obs
echo "chr1	75	200" > exp
check obs exp

echo -e "    general.t40...\c"
$BT merge -i b.bed.gz > obs
echo "chr1	75	200" > exp
check obs exp

echo -e "    general.t41...\c"
$BT merge -i c.bed > obs
echo -n "" > exp
check obs exp

echo -e "    general.t42...\c"
$BT merge -i c.bed.gz > obs
echo -n "" > exp
check obs exp

# allow header lines to start with chr or chrom
echo -e "    general.t43...\c"
echo "chrom	chromStart	chromEnd	name	score	strand
chr14	24800000	24810000	blarg	0	.
chr11	64610000	64620000	blarg	0	.
chr14	24710000	24720000	blarg	0	.
chr11	111230000	111240000	blarg	0	." > exp
$BT intersect -a t.bed -b t.bed -header > obs
check obs exp


rm a.bed.gz b.bed.gz c.bed.gz a.bed b.bed c.bed genome.txt exp obs

[[ $FAILURES -eq 0 ]] || exit 1;
