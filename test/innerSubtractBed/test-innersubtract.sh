ISB=../../bin/innerSubtractBed

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

# cat a.bed
# chr1	100	200	regionA
# chr1	10	400	regionB
# chr1	40	60	regionC
# 
# cat b.bed
# chr1	100	200	regionA
# chr1	10	400	regionB
# chr1	40	60	regionC
# chr1	50	1000	regionA
# chr2	100	500	regionD


###########################################################
# test inner subtraction
###########################################################
echo "    innerSubtractBed.t1...\c"
echo \
"chr1	100	200	regionA
chr1	10	100	regionB
chr1	200	400	regionB" > exp
$ISB -i a.bed --silent > obs
check obs exp
rm obs exp

echo "    innerSubtractBed.t2...\c"
echo \
"chr1	100	200	regionA
chr1	10	100	regionB
chr1	200	400	regionB
chr1	400	1000	regionA
chr2	100	500	regionD" > exp
$ISB -i b.bed --silent > obs
check obs exp
rm obs exp


