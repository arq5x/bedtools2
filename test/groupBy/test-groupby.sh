BT=${BT-../../bin/bedtools}

lines_a=$($BT groupby -g 3-1 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_b=$($BT groupby -g 1-3 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_c=$($BT groupby -g 1,2,3 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_d=$($BT groupby -g 1-2,3 -o collapse -c 4 -i ../map/values3.bed | wc -l)

check(){
    if [ "$1" != "$2" ]; then
        echo "fail groupby" $1 $2
    fi
}

checkfile()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

check $lines_a $lines_b
check $lines_a $lines_c
check $lines_a $lines_d


H=$(head -n 1 values3.header.bed)
A=$($BT groupby -i values3.header.bed -g 1,2,3 -c 4 -o concat -inheader | head -n 1)

if [ "$A" != $'chr1\t0\t10\ta1' ]; then
        echo "fail groupby"
fi

B=$($BT groupby -i values3.header.bed -g 1,2,3 -c 4 -o concat -header | head -n 1)

if [ "$B" != $'#chrom\tstart\tend\tconcat(A)' ]; then
        echo "fail groupby"
fi

###########################################################
#  Test precision
############################################################
echo "    groupby.t01...\c"
echo \
"chr1	11168000	11168003	CALLABLE" > exp
$BT groupby -i test.bed -g 1,4 -c 1,2,3,4 -ops first,first,max,first | cut -f 3-6 > obs
checkfile obs exp

###########################################################
#  Test precision
############################################################
echo "    groupby.t02...\c"
echo \
"chr1	11168000	1.1168e+07	CALLABLE" > exp
$BT groupby -i test.bed -g 1,4 -c 1,2,3,4 -ops first,first,max,first -prec 5 | cut -f 3-6 > obs
checkfile obs exp

rm obs exp