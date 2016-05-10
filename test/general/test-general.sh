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
echo "Error: Invalid record in file -. Record is 
chr1	-1	10" > exp
check obs exp
rm obs exp

###########################################################
#  Enforce start <= end
###########################################################
echo "    general.t02...\c"
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
echo "    general.t03...\c"
echo \
"ERROR: file - has non positional records, which are only valid for the groupBy tool." > exp
echo "chr1	.	2" | $BT merge -i - 2> o
head -n 1 o > obs
check obs exp
rm obs exp



###########################################################
#  Enforce tab-separated files
###########################################################
echo "    general.t05...\c"
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
echo "    general.t06...\c"
$BT merge -i idontexist.bed 2> obs
echo "Error: Unable to open file idontexist.bed. Exiting." > exp
check obs exp
rm obs exp


###########################################################
#  Don't fail on existent, yet empty files.
###########################################################
echo "    general.t07...\c"
$BT merge -i empty.bed 2> obs
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


echo "    general.t09...\c"
$BT intersect -a a.bed -b b.bed > obs
echo "chr1	75	100" > exp
check obs exp

echo "    general.t10...\c"
$BT intersect -a a.bed.gz -b b.bed.gz > obs
echo "chr1	75	100" > exp
check obs exp

echo "    general.t11...\c"
$BT intersect -a a.bed -b b.bed.gz > obs
echo "chr1	75	100" > exp
check obs exp

echo "    general.t12...\c"
$BT intersect -a a.bed.gz -b b.bed > obs
echo "chr1	75	100" > exp
check obs exp

echo "    general.t13...\c"
$BT intersect -a c.bed -b b.bed > obs
echo -n "" > exp
check obs exp

echo "    general.t14...\c"
$BT intersect -a c.bed.gz -b b.bed > obs
echo -n "" > exp
check obs exp

echo "    general.15...\c"
$BT intersect -a c.bed.gz -b c.bed.gz > obs
echo -n "" > exp
check obs exp


echo "    general.t16...\c"
$BT subtract -a a.bed -b b.bed > obs
echo "chr1	1	75" > exp
check obs exp

echo "    general.t17...\c"
$BT subtract -a a.bed.gz -b b.bed.gz > obs
echo "chr1	1	75" > exp
check obs exp

echo "    general.t18...\c"
$BT subtract -a a.bed -b b.bed.gz > obs
echo "chr1	1	75" > exp
check obs exp

echo "    general.t19...\c"
$BT subtract -a a.bed.gz -b b.bed > obs
echo "chr1	1	75" > exp
check obs exp

echo "    general.t20...\c"
$BT subtract -a c.bed -b b.bed > obs
echo -n "" > exp
check obs exp

echo "    general.t21...\c"
$BT subtract -a c.bed.gz -b b.bed > obs
echo -n "" > exp
check obs exp

echo "    general.22...\c"
$BT subtract -a c.bed.gz -b c.bed.gz > obs
echo -n "" > exp
check obs exp


echo "    general.t23...\c"
$BT window -a a.bed -b b.bed > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo "    general.t24...\c"
$BT window -a a.bed.gz -b b.bed.gz > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo "    general.t25...\c"
$BT window -a a.bed -b b.bed.gz > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo "    general.t26...\c"
$BT window -a a.bed.gz -b b.bed > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo "    general.t27...\c"
$BT window -a c.bed -b b.bed > obs
echo -n "" > exp
check obs exp

echo "    general.t28...\c"
$BT window -a c.bed.gz -b b.bed > obs
echo -n "" > exp
check obs exp

echo "    general.29...\c"
$BT window -a c.bed.gz -b c.bed.gz > obs
echo -n "" > exp
check obs exp


echo "    general.t30...\c"
$BT closest -a a.bed -b b.bed > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo "    general.t31...\c"
$BT closest -a a.bed.gz -b b.bed.gz > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo "    general.t32...\c"
$BT closest -a a.bed -b b.bed.gz > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo "    general.t33...\c"
$BT closest -a a.bed.gz -b b.bed > obs
echo "chr1	1	100	chr1	75	200" > exp
check obs exp

echo "    general.t34...\c"
$BT closest -a c.bed -b b.bed > obs
echo -n "" > exp
check obs exp

echo "    general.t35...\c"
$BT closest -a c.bed.gz -b b.bed > obs
echo -n "" > exp
check obs exp

echo "    general.36...\c"
$BT closest -a c.bed.gz -b c.bed.gz > obs
echo -n "" > exp
check obs exp


echo "    general.t37...\c"
$BT merge -i a.bed > obs
echo "chr1	1	100" > exp
check obs exp

echo "    general.t38...\c"
$BT merge -i a.bed.gz > obs
echo "chr1	1	100" > exp
check obs exp

echo "    general.t39...\c"
$BT merge -i b.bed > obs
echo "chr1	75	200" > exp
check obs exp

echo "    general.t40...\c"
$BT merge -i b.bed.gz > obs
echo "chr1	75	200" > exp
check obs exp

echo "    general.t41...\c"
$BT merge -i c.bed > obs
echo -n "" > exp
check obs exp

echo "    general.t42...\c"
$BT merge -i c.bed.gz > obs
echo -n "" > exp
check obs exp

rm a.bed.gz b.bed.gz c.bed.gz a.bed b.bed c.bed genome.txt exp obs

