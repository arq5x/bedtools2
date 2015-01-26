BT=${BT-../../bin/bedtools}


cat /dev/zero |tr "\0" "\n" | head -n 10000 |\
awk '{S=int(rand()*1000.0);E=S+int(rand()*100000); printf("chrX\t%d\t%d\n",S,E);}' > tmp.bed

echo "    slpit.size...\c"
rm -f _tmp.*.bed
${BT} split -i tmp.bed -p _tmp -n 50 -a size > _tmp.size.tsv
cat _tmp.size.tsv


echo "    slpit.simple...\c"
rm -f _tmp.*.bed
${BT} split -i tmp.bed -p _tmp -n 50 -a simple > _tmp.simple.tsv
cat _tmp.simple.tsv

echo "    slpit.simple...\c"
rm -f _tmp.*.bed
${BT} split -i tmp.bed -p _tmp -n 1000 -a simple > _tmp.simple.tsv
cat _tmp.simple.tsv



rm -f _tmp.*.bed jeter.bed tmp.bed  _tmp.simple.tsv  _tmp.size.tsv
