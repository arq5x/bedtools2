echo " Performing general tests:"
cd general; sh test-general.sh; cd ..

echo " Testing bedtools bed12tobed6:"
cd bed12tobed6; sh test-bed12tobed6.sh; cd ..

echo " Testing bedtools bamtobed:"
cd bamtobed; sh test-bamtobed.sh; cd ..

echo " Testing bedtools closest:"
cd closest; sh test-closest.sh; cd ..

echo " Testing bedtools cluster:"
cd cluster; sh test-cluster.sh; cd ..

echo " Testing bedtools coverage:"
cd coverage; sh test-coverage.sh; cd ..

echo " Testing bedtools expand:"
cd expand; sh test-expand.sh; cd ..

echo " Testing bedtools flank:"
cd flank; sh test-flank.sh; cd ..

echo " Testing bedtools genomecov:"
cd genomecov; sh test-genomecov.sh; cd ..

echo " Testing bedtools intersect:"
cd intersect; sh test-intersect.sh; cd ..

echo " Testing bedtools map:"
cd map; sh test-map.sh; cd ..

echo " Testing bedtools merge:"
cd merge; sh test-merge.sh; cd ..

echo " Testing bedtools slop:"
cd slop; sh test-slop.sh; cd ..
