echo " Performing general tests:"
cd general; bash test-general.sh; cd ..

echo " Testing bedtools bed12tobed6:"
cd bed12tobed6; bash test-bed12tobed6.sh; cd ..

echo " Testing bedtools bamtobed:"
cd bamtobed; bash test-bamtobed.sh; cd ..

echo " Testing bedtools closest:"
cd closest; bash test-closest.sh; cd ..

echo " Testing bedtools cluster:"
cd cluster; bash test-cluster.sh; cd ..

echo " Testing bedtools coverage:"
cd coverage; bash test-coverage.sh; cd ..

echo " Testing bedtools expand:"
cd expand; bash test-expand.sh; cd ..

echo " Testing bedtools flank:"
cd flank; bash test-flank.sh; cd ..

echo " Testing bedtools fisher:"
cd fisher; bash test-fisher.sh; cd ..

echo " Testing bedtools genomecov:"
cd genomecov; bash test-genomecov.sh; cd ..

echo " Testing bedtools getfasta:"
cd getfasta; bash test-getfasta.sh; cd ..

echo " Testing bedtools intersect:"
cd intersect; bash test-intersect.sh; bash new_test-intersect.sh; cd ..

echo " Testing bedtools jaccard:"
cd jaccard; bash test-jaccard.sh; cd ..

echo " Testing bedtools map:"
cd map; bash test-map.sh; cd ..

echo " Testing bedtools merge:"
cd merge; bash test-merge.sh; cd ..

echo " Testing bedtools multicov:"
cd multicov; bash test-multicov.sh; cd ..

echo " Testing bedtools reldist:"
cd reldist; bash test-reldist.sh; cd ..

echo " Testing bedtools shift:"
cd shift; bash test-shift.sh; cd ..

echo " Testing bedtools slop:"
cd slop; bash test-slop.sh; cd ..

echo " Testing bedtools sort:"
cd slop; bash test-sort.sh; cd ..

echo " Testing bedtools shuffle:"
cd shuffle; bash test-shuffle.sh; cd ..

echo " Testing bedtools subtract:"
cd subtract; bash test-subtract.sh; cd ..

echo " Testing bedtools sample:"
cd sample; bash test-sample.sh; cd ..

echo " Testing bedtools split:"
cd split; bash test-split.sh; cd ..

echo " Testing bedtools spacing:"
cd spacing; bash test-spacing.sh; cd ..
