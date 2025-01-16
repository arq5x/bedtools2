set -e;
BT=${BT-../../bin/bedtools}

FAILURES=0;

check()
{
     if diff -Z $1 $2; then
          echo ok
     else
          FAILURES=$(expr $FAILURES + 1);
          echo fail
     fi
}


###########################################################
# test baseline complement
###########################################################
echo -e "    complement.t1...\c"
echo "chr1	1	20" > exp
echo -e "chr1\t0\t1" > tmp1
echo -e "chr1\t20"   > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               > obs
check obs exp
rm tmp1 tmp2
rm obs exp


###########################################################
# ends are covered
###########################################################
echo -e "    complement.t2...\c"
echo "chr1	1	19" > exp
echo -e "chr1\t0\t1\nchr1\t19\t20" > tmp1
echo -e "chr1\t20"                 > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               > obs
check obs exp
rm tmp1 tmp2
rm obs exp

###########################################################
# middle is covered
###########################################################
echo -e "    complement.t3...\c"
echo "chr1	0	10
chr1	15	20" > exp
echo -e "chr1\t10\t15" > tmp1
echo -e "chr1\t20"     > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               > obs
check obs exp
rm tmp1 tmp2
rm obs exp

###########################################################
# entirety is covered
###########################################################
echo -e "    complement.t4...\c"
touch exp
echo -e "chr1\t0\t20" > tmp1
echo -e "chr1\t20"    > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               > obs
check obs exp
rm tmp1 tmp2
rm obs exp


###########################################################
# nothing is covered
###########################################################
echo -e "    complement.t5...\c"
echo "chr1	0	20" > exp
echo -e "chr2\t0\t20"          > tmp1
echo -e "chr1\t20\nchr2\t20"   > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               > obs
check obs exp
rm tmp1 tmp2
rm obs exp

###########################################################
# Issue #356
###########################################################
echo -e "    complement.t6...\c"
echo "chr1	10000	249250621" > exp
echo -e "chr1\t0\t10000\ttelomere" > tmp1
echo -e "chr1\t249250621"          > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               > obs
check obs exp
rm tmp1 tmp2
rm obs exp

###########################################################
# Multiple chroms 
###########################################################
echo -e "    complement.t7...\c"
echo "chr1	0	10
chr2	0	10" > exp
echo -e "chr1\t10\t20\nchr2\t10\t20" > tmp1
echo -e "chr1\t20\nchr2\t20"         > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               > obs
check obs exp
rm tmp1 tmp2
rm obs exp

###########################################################
# Multiple chroms, chr1 is covered
###########################################################
echo -e "    complement.t8...\c"
echo "chr2	0	10" > exp
echo -e "chr1\t0\t20\nchr2\t10\t20" > tmp1
echo -e "chr1\t20\nchr2\t20"        > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               > obs
check obs exp
rm tmp1 tmp2
rm obs exp

###########################################################
# record exceeds chrom length
###########################################################
echo -e "    complement.t9...\c"
echo -e "***** WARNING: chr1:90-110 exceeds the length of chromosome (chr1)\nchr1\t0\t90" > exp
echo -e "chr1\t90\t110" > tmp1
echo -e "chr1\t100"     > tmp2
$BT complement -i tmp1 \
               -g tmp2 \
               &> obs
check obs exp
rm tmp1 tmp2
rm obs exp
[[ $FAILURES -eq 0 ]] || exit 1;

###########################################################
# -L only reports chroms that were in the BED file.
###########################################################
echo -e "    complement.t9...\c"
echo "chr1	0	1
chr1	500	900
chr1	950	1000" > exp
$BT complement -i issue_503.bed \
               -g issue_503.genome \
               -L \
               > obs
check obs exp
rm obs exp

###########################################################
# Now, without -L
###########################################################
echo -e "    complement.t10...\c"
echo "chr1	0	1
chr1	500	900
chr1	950	1000
chr2	0	1000
chr3	0	1000" > exp
$BT complement -i issue_503.bed \
               -g issue_503.genome \
               > obs
check obs exp
rm obs exp

[[ $FAILURES -eq 0 ]] || exit 1;
