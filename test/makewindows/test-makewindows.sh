BT=${BT-../../bin/bedtools}

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}


###########################################################
#  Test window size / forward window numbering 
############################################################
echo "    makewindows.t01...\c"
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
echo "    makewindows.t02...\c"
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
echo "    makewindows.t03...\c"
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
echo "    makewindows.t04...\c"
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
echo "    makewindows.t05...\c"
echo \
"chr5	60000	63334	1
chr5	63334	66668	2
chr5	66668	70000	3
chr5	73000	78667	1
chr5	78667	84334	2
chr5	84334	90000	3
chr5	100000	100334	1
chr5	100334	100668	2
chr5	100668	101000	3" > exp
$BT makewindows -n 3 -b input.bed -i winnum > obs
check obs exp
rm obs exp

###########################################################
#  Test fixed size / reverse window numbering 
############################################################
echo "    makewindows.t06...\c"
echo \
"chr5	60000	63334	3
chr5	63334	66668	2
chr5	66668	70000	1
chr5	73000	78667	3
chr5	78667	84334	2
chr5	84334	90000	1
chr5	100000	100334	3
chr5	100334	100668	2
chr5	100668	101000	1" > exp
$BT makewindows -n 3 -b input.bed -reverse -i winnum > obs
check obs exp
rm obs exp
