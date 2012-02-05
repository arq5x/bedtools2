echo " Performing general tests:"
cd general
sh test-general.sh
cd ..


echo " Testing bedtools intersect:"
cd intersect
sh test-intersect.sh
cd ..

echo " Testing bedtools map:"
cd map
sh test-map.sh
cd ..

echo " Testing bedtools merge:"
cd merge
sh test-merge.sh
cd ..