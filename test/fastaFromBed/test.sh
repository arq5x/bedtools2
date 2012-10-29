LINES=$(echo $'chr1\t1\t10' | ../../bin/bedtools getfasta -fi t.fa -bed stdin -fo - | awk 'END{ print NR }')


if [ "$LINES" != "2" ]; then
    echo "BAD";
fi

LEN=$(../../bin/bedtools getfasta -split -fi t.fa -bed blocks.bed -fo stdout | awk '(NR == 2){ print length($0) }')

if [ "$LEN" != "22" ]; then
    echo "BAD";
fi

SEQ=$(../../bin/bedtools getfasta -split -fi t.fa -bed blocks.bed -fo stdout | awk '(NR == 4){ print $0 }')

if [ "$SEQ" != "cta" ]; then
    echo "BAD";
fi

