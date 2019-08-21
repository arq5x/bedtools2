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

# Tools that are expected to fail should exit non-zero. If they don't, that is a
# failure.

failedtofail()
{
	FAILURES=$(expr $FAILURES + 1);
	echo "Expected non-zero exit status but didn't get one";
}

###########################################################
#  Enforce non-negative coordinates
###########################################################
echo -e "    general.t01...\c"
echo \
"chr1	1	10
chr1	-1	10" | $BT merge -i - 2> obs \
	&& failedtofail || true;
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
chr1	10	5" | $BT merge -i - 2> obs \
	&& failedtofail || true;
echo "Error: Invalid record in file -. Record is 
chr1	10	5" > exp
check obs exp
rm obs exp

###########################################################
#  Enforce integer coordinates
###########################################################
echo -e "    general.t03...\c"

echo "Error: unable to open file or unable to determine types for file -

- Please ensure that your file is TAB delimited (e.g., cat -t FILE).
- Also ensure that your file has integer chromosome coordinates in the 
  expected columns (e.g., cols 2 and 3 for BED)." >> exp

#echo "ERROR: file - has non positional records, which are only valid for " >> exp
#echo " the groupBy tool. Perhaps you are using a header line(s) that starts with " >> exp
#echo " something other than \"#\", \"chrom\", or \"chr\" (any case)?" >> exp

echo "chr1	.	2" | $BT merge -i - 2> o \
	&& failedtofail || true;
tail -n 5 o > obs
check obs exp
rm o obs exp



###########################################################
#  Enforce tab-separated files
###########################################################
echo -e "    general.t05...\c"
echo \
"chr1 1 2" | $BT merge -i - 2> o \
	&& failedtofail || true;
echo "Error: unable to open file or unable to determine types for file -

- Please ensure that your file is TAB delimited (e.g., cat -t FILE).
- Also ensure that your file has integer chromosome coordinates in the 
  expected columns (e.g., cols 2 and 3 for BED)." >> exp
tail -n 5 o > obs
check obs exp
rm o obs exp

###########################################################
#  Fail on non-existent files.
###########################################################
echo -e "    general.t06...\c"
$BT merge -i idontexist.bed 2> obs \
	&& failedtofail || true;
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

echo -e "    general.t44...\c"
$BT intersect -a a.bed.gz -b b.bed -g hg19.fa.fai > obs
echo "chr1	75	100" > exp
check obs exp

# allow header lines to start with chrom or track
echo -e "    general.t45...\c"
echo "chr1	15	20" > exp
$BT intersect -a a.trackheader.bed -b a.chromheader.bed > obs
check obs exp

# allow header lines to start with chrom or track
echo -e "    general.t46...\c"
echo "track
chr1	15	20" > exp
$BT intersect -a a.trackheader.bed -b a.chromheader.bed -header > obs
check obs exp


rm a.bed.gz b.bed.gz c.bed.gz a.bed b.bed c.bed genome.txt exp obs

[[ $FAILURES -eq 0 ]] || exit 1;
