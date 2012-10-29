# t.fa
# >chr1
# aggggggggg
# cggggggggg
# tggggggggg
# aggggggggg
# cggggggggg

check()
{
	if $1 = $2; then
    	echo ok
	else
    	echo fail
	fi
}

echo "    getfasta.t01...\c"
LINES=$(echo $'chr1\t1\t10' | ../../bin/bedtools getfasta -fi t.fa -bed stdin -fo - | awk 'END{ print NR }')
if [ "$LINES" != "2" ]; then
    echo fail
else
    echo ok
fi

echo "    getfasta.t02...\c"
LEN=$(../../bin/bedtools getfasta -split -fi t.fa -bed blocks.bed -fo stdout | awk '(NR == 2){ print length($0) }')
if [ "$LINES" != "2" ]; then
    echo fail
else
    echo ok
fi

echo "    getfasta.t03...\c"
SEQ=$(../../bin/bedtools getfasta -split -fi t.fa -bed blocks.bed -fo stdout | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "cta" ]; then
    echo fail
else
    echo ok
fi

# test -fo -
echo "    getfasta.t04...\c"
SEQ=$(../../bin/bedtools getfasta -split -fi t.fa -bed blocks.bed -fo - | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "cta" ]; then
    echo fail
else
    echo ok
fi


# test -split with -s -
echo "    getfasta.t05...\c"
SEQ=$(../../bin/bedtools getfasta -split -s -fi t.fa -bed blocks.bed -fo - | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "tag" ]; then
    echo fail
else
    echo ok
fi