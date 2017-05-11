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
#  Test window size / forward window numbering 
############################################################
echo -e "    makewindows.t01...\c"
echo \
"chr5	60000	65000	1
chr5	65000	70000	2
chr5	73000	78000	1
chr5	78000	83000	2
chr5	83000	88000	3
chr5	88000	90000	4
chr5	100000	101000	1" > exp
$BT makewindows -b input.bed -w 5000  -i winnum > obs
check obs exp
rm obs exp

###########################################################
#  Test window size / reverse window numbering 
############################################################
echo -e "    makewindows.t02...\c"
echo \
"chr5	60000	65000	2
chr5	65000	70000	1
chr5	73000	78000	4
chr5	78000	83000	3
chr5	83000	88000	2
chr5	88000	90000	1
chr5	100000	101000	1" > exp
$BT makewindows -b input.bed -reverse -w 5000 -i winnum > obs
check obs exp
rm obs exp

###########################################################
#  Test window+step size / forward window numbering 
############################################################
echo -e "    makewindows.t03...\c"
echo \
"chr5	60000	65000	1
chr5	62000	67000	2
chr5	64000	69000	3
chr5	66000	70000	4
chr5	68000	70000	5
chr5	73000	78000	1
chr5	75000	80000	2
chr5	77000	82000	3
chr5	79000	84000	4
chr5	81000	86000	5
chr5	83000	88000	6
chr5	85000	90000	7
chr5	87000	90000	8
chr5	89000	90000	9
chr5	100000	101000	1" > exp
$BT makewindows -b input.bed -w 5000 -s 2000 -i winnum > obs
check obs exp
rm obs exp

###########################################################
#  Test window size / reverse window numbering 
############################################################
echo -e "    makewindows.t04...\c"
echo \
"chr5	60000	65000	5
chr5	62000	67000	4
chr5	64000	69000	3
chr5	66000	70000	2
chr5	68000	70000	1
chr5	73000	78000	9
chr5	75000	80000	8
chr5	77000	82000	7
chr5	79000	84000	6
chr5	81000	86000	5
chr5	83000	88000	4
chr5	85000	90000	3
chr5	87000	90000	2
chr5	89000	90000	1
chr5	100000	101000	1" > exp
$BT makewindows -b input.bed -reverse -w 5000 -s 2000 -i winnum > obs
check obs exp
rm obs exp


###########################################################
#  Test fixed size / forward window numbering 
############################################################
echo -e "    makewindows.t05...\c"
echo \
"chr5	60000	63333	1
chr5	63333	66666	2
chr5	66666	70000	3
chr5	73000	78666	1
chr5	78666	84332	2
chr5	84332	90000	3
chr5	100000	100333	1
chr5	100333	100666	2
chr5	100666	101000	3" > exp
$BT makewindows -n 3 -b input.bed -i winnum > obs
check obs exp
rm obs exp

###########################################################
#  Test fixed size / reverse window numbering 
############################################################
echo -e "    makewindows.t06...\c"
echo \
"chr5	60000	63333	3
chr5	63333	66666	2
chr5	66666	70000	1
chr5	73000	78666	3
chr5	78666	84332	2
chr5	84332	90000	1
chr5	100000	100333	3
chr5	100333	100666	2
chr5	100666	101000	1" > exp
$BT makewindows -n 3 -b input.bed -reverse -i winnum > obs
check obs exp
rm obs exp

###########################################################
#  Test that we alway get the number of requested windows
###########################################################
echo -e "    makewindows.t07...\c"
echo \
"1	11	14	A_1
1	14	17	A_2
1	17	20	A_3
1	20	23	A_4
1	23	26	A_5
1	26	29	A_6
1	29	32	A_7
1	32	35	A_8
1	35	38	A_9
1	38	44	A_10" > exp
$BT makewindows -b a.33bp.bed -n 10 -i srcwinnum> obs
check obs exp
rm obs exp

###########################################################
#  Test that we alway get the number of requested windows
###########################################################
echo -e "    makewindows.t08...\c"
echo \
"1	10	12	B_1
1	12	14	B_2
1	14	16	B_3
1	16	18	B_4
1	18	20	B_5
1	20	22	B_6
1	22	24	B_7
1	24	26	B_8
1	26	28	B_9
1	28	30	B_10" > exp
$BT makewindows -b a.20bp.bed -n 10 -i srcwinnum> obs
check obs exp
rm obs exp

###########################################################
#  interval is smaller than n
###########################################################
echo -e "    makewindows.t09...\c"
echo \
"WARNING: Interval 1:10-19 is smaller than the number of windows requested. Skipping." > exp
$BT makewindows -b a.19bp.bed -n 10 -i srcwinnum 2> obs
check obs exp
rm obs exp
[[ $FAILURES -eq 0 ]] || exit 1;
