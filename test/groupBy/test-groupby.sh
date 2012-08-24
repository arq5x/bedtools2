lines_a=$(../../bin/groupBy -g 3-1 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_b=$(../../bin/groupBy -g 1-3 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_c=$(../../bin/groupBy -g 1,2,3 -o collapse -c 4 -i ../map/values3.bed | wc -l)
lines_d=$(../../bin/groupBy -g 1-2,3 -o collapse -c 4 -i ../map/values3.bed | wc -l)

check(){
    if [ "$1" != "$2" ]; then
        "fail groupby"
    fi
}

check $lines_a $lines_b
check $lines_a $lines_c
check $lines_a $lines_d
