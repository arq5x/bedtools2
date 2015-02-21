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
#  Test forward window numbering 
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
#  Test reverse window numbering 
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
$BT makewindows -b input.bed -reverse -w 5000  -i winnum > obs

check obs exp
rm obs exp
