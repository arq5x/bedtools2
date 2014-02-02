# t.fa
# >chr1
# aggggggggg
# cggggggggg
# tggggggggg
# aggggggggg
# cggggggggg

BT=${BT-../../bin/bedtools}

check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
    fi
}

echo "    getfasta.t01...\c"
LINES=$(echo $'chr1\t1\t10' | $BT getfasta -fi t.fa -bed stdin -fo - | awk 'END{ print NR }')
if [ "$LINES" != "2" ]; then
    echo fail
else
    echo ok
fi

echo "    getfasta.t02...\c"
LEN=$($BT getfasta -split -fi t.fa -bed blocks.bed -fo stdout | awk '(NR == 2){ print length($0) }')
if [ "$LINES" != "2" ]; then
    echo fail
else
    echo ok
fi

echo "    getfasta.t03...\c"
SEQ=$($BT getfasta -split -fi t.fa -bed blocks.bed -fo stdout | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "cta" ]; then
    echo fail
else
    echo ok
fi

# test -fo -
echo "    getfasta.t04...\c"
SEQ=$($BT getfasta -split -fi t.fa -bed blocks.bed -fo - | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "cta" ]; then
    echo fail
else
    echo ok
fi


# test -split with -s -
echo "    getfasta.t05...\c"
SEQ=$($BT getfasta -split -s -fi t.fa -bed blocks.bed -fo - | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "tag" ]; then
    echo fail
else
    echo ok
fi

# test -fullHeader
echo "    getfasta.t06...\c"
LINES=$(echo $'chr1 assembled by consortium X\t1\t10' | $BT getfasta -fullHeader -fi t_fH.fa -bed stdin -fo - | awk 'END{ print NR }')
if [ "$LINES" != "2" ]; then
    echo fail
else
    echo ok
fi

# test without -fullHeader
echo "    getfasta.t07...\c"
echo "WARNING. chromosome (chr1 assembled by consortium X) was not found in the FASTA file. Skipping." > exp
echo $'chr1 assembled by consortium X\t1\t10' | $BT getfasta -fi t_fH.fa -bed - -fo - 2> obs

check obs exp

rm obs exp