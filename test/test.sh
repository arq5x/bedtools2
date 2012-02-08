echo " Performing general tests:"
cd general
sh test-general.sh
cd ..

echo " Testing bedtools bed12tobed6:"
cd bed12tobed6
sh test-bed12tobed6.sh
cd ..

echo " Testing bedtools bamtobed:"
cd bamtobed
sh test-bamtobed.sh
cd ..

echo " Testing bedtools genomecov:"
cd genomecov
sh test-genomecov.sh
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