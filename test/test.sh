echo " Testing bedtools intersect:"
cd intersect
sh test-intersect.sh
cd ..

echo " Testing bedtools merge:"
cd merge
sh test-merge.sh
cd ..