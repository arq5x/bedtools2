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
#  Test forward and reverse window numbering 
############################################################
echo "    makewindows.t01...\c"
echo

$BT makewindows -b input.bed -w 500 -i winnum > obs
check obs exp
rm obs exp
