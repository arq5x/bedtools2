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
# test baseline complement
###########################################################
echo "    complement.t1...\c"
echo "chr1	1	20" > exp
$BT complement -i <(echo -e "chr1\t0\t1") \
               -g <(echo -e "chr1\t20") \
               > obs
check obs exp
rm obs exp


###########################################################
# ends are covered
###########################################################
echo "    complement.t2...\c"
echo "chr1	1	19" > exp
$BT complement -i <(echo -e "chr1\t0\t1\nchr1\t19\t20") \
               -g <(echo -e "chr1\t20") \
               > obs
check obs exp
rm obs exp

###########################################################
# middle is covered
###########################################################
echo "    complement.t3...\c"
echo "chr1	0	10
chr1	15	20" > exp
$BT complement -i <(echo -e "chr1\t10\t15") \
               -g <(echo -e "chr1\t20") \
               > obs
check obs exp
rm obs exp

###########################################################
# entirety is covered
###########################################################
echo "    complement.t4...\c"
touch exp
$BT complement -i <(echo -e "chr1\t0\t20") \
               -g <(echo -e "chr1\t20") \
               > obs
check obs exp
rm obs exp


###########################################################
# nothing is covered
###########################################################
echo "    complement.t5...\c"
echo "chr1	0	20" > exp
$BT complement -i <(echo -e "chr2\t0\t20") \
               -g <(echo -e "chr1\t20\nchr2\t20") \
               > obs
check obs exp
rm obs exp

###########################################################
# Issue #356
###########################################################
echo "    complement.t6...\c"
echo "chr1	10000	249250621" > exp
$BT complement -i <(echo -e "chr1\t0\t10000\ttelomere") \
               -g <(echo -e "chr1\t249250621") \
               > obs
check obs exp
rm obs exp

###########################################################
# Multiple chroms 
###########################################################
echo "    complement.t7...\c"
echo "chr1	0	10
chr2	0	10" > exp
$BT complement -i <(echo -e "chr1\t10\t20\nchr2\t10\t20") \
               -g <(echo -e "chr1\t20\nchr2\t20") \
               > obs
check obs exp
rm obs exp

###########################################################
# Multiple chroms, chr1 is covered
###########################################################
echo "    complement.t8...\c"
echo "chr2	0	10" > exp
$BT complement -i <(echo -e "chr1\t0\t20\nchr2\t10\t20") \
               -g <(echo -e "chr1\t20\nchr2\t20") \
               > obs
check obs exp
rm obs exp

###########################################################
# record exceeds chrom length
###########################################################
echo "    complement.t9...\c"
echo -e "***** WARNING: chr1:90-110 exceeds the length of chromosome (chr1)\nchr1\t0\t90" > exp
$BT complement -i <(echo -e "chr1\t90\t110") \
               -g <(echo -e "chr1\t100") \
               &> obs
check obs exp
rm obs exp
