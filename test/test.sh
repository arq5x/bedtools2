# In the event of a test failing, exit with non-zero status:
set -e;

STARTWD=$(pwd);
FAILURES=0;

echo " Performing general tests:"
(cd $STARTWD/general && bash test-general.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools bed12tobed6:"
(cd $STARTWD/bed12tobed6 && bash test-bed12tobed6.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools bamtobed:"
(cd $STARTWD/bamtobed && bash test-bamtobed.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools closest:"
(cd $STARTWD/closest && bash test-closest.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools cluster:"
(cd $STARTWD/cluster && bash test-cluster.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools coverage:"
(cd $STARTWD/coverage && bash test-coverage.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools expand:"
(cd $STARTWD/expand && bash test-expand.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools flank:"
(cd $STARTWD/flank && bash test-flank.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools fisher:"
(cd $STARTWD/fisher && bash test-fisher.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools genomecov:"
(cd $STARTWD/genomecov && bash test-genomecov.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools getfasta:"
(cd $STARTWD/getfasta && bash test-getfasta.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools groupby:"
(cd $STARTWD/groupby && bash test-groupby.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools intersect:"
(cd $STARTWD/intersect && bash test-intersect.sh; bash new_test-intersect.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools jaccard:"
(cd $STARTWD/jaccard && bash test-jaccard.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools map:"
(cd $STARTWD/map && bash test-map.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools merge:"
(cd $STARTWD/merge && bash test-merge.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools multicov:"
(cd $STARTWD/multicov && bash test-multicov.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools reldist:"
(cd $STARTWD/reldist && bash test-reldist.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools shift:"
(cd $STARTWD/shift && bash test-shift.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools slop:"
(cd $STARTWD/slop && bash test-slop.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools sort:"
(cd $STARTWD/sort && bash test-sort.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools shuffle:"
(cd $STARTWD/shuffle && bash test-shuffle.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools subtract:"
(cd $STARTWD/subtract && bash test-subtract.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools sample:"
(cd $STARTWD/sample && bash test-sample.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools split:"
(cd $STARTWD/split && bash test-split.sh) || FAILURES=$(expr $FAILURES + 1);

echo " Testing bedtools spacing:"
(cd $STARTWD/spacing && bash test-spacing.sh) || FAILURES=$(expr $FAILURES + 1);

if [[ $FAILURES -gt 0 ]]; then
    >&2 echo "$FAILURES tests failed";
    exit 1;
fi



