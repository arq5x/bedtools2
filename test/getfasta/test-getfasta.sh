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
LINES=$(echo $'chr1\t1\t10' | $BT getfasta -fi t.fa -bed stdin | awk 'END{ print NR }')
if [ "$LINES" != "2" ]; then
    echo fail
else
    echo ok
fi

echo "    getfasta.t02...\c"
LEN=$($BT getfasta -split -fi t.fa -bed blocks.bed | awk '(NR == 2){ print length($0) }')
if [ "$LINES" != "2" ]; then
    echo fail
else
    echo ok
fi

echo "    getfasta.t03...\c"
SEQ=$($BT getfasta -split -fi t.fa -bed blocks.bed | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "cta" ]; then
    echo fail
else
    echo ok
fi

# test -fo -
echo "    getfasta.t04...\c"
SEQ=$($BT getfasta -split -fi t.fa -bed blocks.bed | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "cta" ]; then
    echo fail
else
    echo ok
fi


# test -split with -s -
echo "    getfasta.t05...\c"
SEQ=$($BT getfasta -split -s -fi t.fa -bed blocks.bed | awk '(NR == 4){ print $0 }')
if [ "$SEQ" != "tag" ]; then
    echo fail
else
    echo ok
fi

# test -fullHeader
echo "    getfasta.t06...\c"
LINES=$(echo $'chr1 assembled by consortium X\t1\t10' | $BT getfasta -fullHeader -fi t_fH.fa -bed stdin | awk 'END{ print NR }')
if [ "$LINES" != "2" ]; then
    echo fail
else
    echo ok
fi

# test without -fullHeader
echo "    getfasta.t07...\c"
echo "WARNING. chromosome (chr1 assembled by consortium X) was not found in the FASTA file. Skipping." > exp
echo $'chr1 assembled by consortium X\t1\t10' | $BT getfasta -fi t_fH.fa -bed - 2> obs

check obs exp
rm obs exp


# test IUPAC
echo "    getfasta.t08...\c"
echo \
">1:0-16
AGCTYRWSKMDVHBXN
>2:0-16
agctyrwskmdvhbxn" > exp
$BT getfasta  -fi test.iupac.fa -bed test.iupac.bed > obs
check obs exp
rm obs exp test.iupac.fa.fai


# test IUPAC revcomp
echo "    getfasta.t09...\c"
echo \
">1:0-16(-)
NXVDBHKMSWYRAGCT
>2:0-16(-)
nxvdbhkmswyragct" > exp
$BT getfasta  -fi test.iupac.fa -bed test.iupac.bed -s > obs
check obs exp
rm obs exp test.iupac.fa.fai

# test the warning about an outdated FASTA index file
echo "    getfasta.t10...\c"
echo \
">chr1
cggggggggg
>chr2
AAATTTTTTTTTT" > test.fa
# create an index file
echo -e "chr2\t2\t10" | $BT getfasta -fi test.fa -bed - > /dev/null
# modify the FASTA file in a second
sleep 1 
touch test.fa
echo -e "chr2\t2\t10" | $BT getfasta -fi test.fa -bed -  \
	> /dev/null 2> obs
echo "Warning: the index file is older than the FASTA file." > exp
check obs exp
rm obs exp test.fa test.fa.fai


echo "    getfasta.t11...\c"
echo ">chr1:0-40
agggggggggcgggggggggtgggggggggaggggggggg
>chr1:0-40
agggggggggcgggggggggtgggggggggaggggggggg" > exp
$BT getfasta -fi t.fa -bed blocks.bed > obs
check obs exp
rm obs exp

echo "    getfasta.t12...\c"
echo ">three_blocks_match::chr1:0-40
agggggggggcgggggggggtgggggggggaggggggggg
>three_blocks_match::chr1:0-40
agggggggggcgggggggggtgggggggggaggggggggg" > exp
$BT getfasta -fi t.fa -bed blocks.bed -name > obs
check obs exp
rm obs exp

echo "    getfasta.t13...\c"
echo ">three_blocks_match::chr1:0-40(+)
agggggggggcgggggggggtgggggggggaggggggggg
>three_blocks_match::chr1:0-40(-)
ccccccccctcccccccccacccccccccgccccccccct" > exp
$BT getfasta -fi t.fa -bed blocks.bed -name -s > obs
check obs exp
rm obs exp


echo "    getfasta.t134..\c"
echo ">chr1:0-1
a" > exp
$BT getfasta -fi t.fa -bed <(echo -e "chr1\t0\t1") > obs
check obs exp
rm obs exp