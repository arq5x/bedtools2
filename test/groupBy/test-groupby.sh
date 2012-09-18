lines_a=$(../../bin/groupBy -g 3-1 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_b=$(../../bin/groupBy -g 1-3 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_c=$(../../bin/groupBy -g 1,2,3 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_d=$(../../bin/groupBy -g 1-2,3 -o collapse -c 4 -i ../map/values3.bed | wc -l)

check(){
    if [ "$1" != "$2" ]; then
        echo "fail groupby" $1 $2
    fi
}

check $lines_a $lines_b
check $lines_a $lines_c
check $lines_a $lines_d


H=$(head -n 1 values3.header.bed)
A=$(../../bin/bedtools groupby -i values3.header.bed -g 1,2,3 -c 4 -o concat -inheader | head -n 1)

if [ "$A" != $'chr1\t0\t10\ta1' ]; then
        echo "fail groupby"
fi

B=$(../../bin/bedtools groupby -i values3.header.bed -g 1,2,3 -c 4 -o concat -header | head -n 1)

if [ "$B" != $'#chrom\tstart\tend\tconcat(A)' ]; then
        echo "fail groupby"
fi
