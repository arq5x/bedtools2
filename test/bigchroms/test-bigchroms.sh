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
#  Test a basic self intersection
############################################################
echo -e "    bigchroms.t01...\c"
$BT intersect -sorted -a abig.bed -b abig.bed > obs
check obs abig.bed
rm obs

echo -e "    bigchroms.t02...\c"
$BT intersect -a abig.bed -b abig.bed > obs
check obs abig.bed
rm obs

if [[ "$BT_NO_BIG_FILES" != "" ]]; then
python make-big-chrom.py

echo -e "    bigchroms.t03...big get fasta \c"
$BT getfasta -fi bigx.fasta -bed bigx.bed | tail -1 > obs
echo "ACTGACCCCGAGACGTTTGCATCCTGCACAGCTAGAGATCCTTTATTAAAAGCACACTGT" > exp
check obs exp
rm obs exp

rm bigx.fasta*

else

echo -e "    bigchroms.t03... skip \c"
echo " set env var 'BT_NO_BIG_FILES' to run this test"
fi

echo -e "    bigchroms.t04... merge \c"

echo "chr1	1	9000000000" > exp
$BT merge -i big4.bed > obs
check obs exp
rm obs exp


echo -e "    bigchroms.t05... closest \c"
$BT closest -a big4c.bed -b big4.bed  > obs

echo "chr1	1	10	chr1	1	10
chr1	1	10	chr1	9	8000000000
chr1	8900000000	9000000000	chr1	8000000000	9000000000" > exp
check obs exp
rm obs exp

[[ $FAILURES -eq 0 ]] || exit 1;


