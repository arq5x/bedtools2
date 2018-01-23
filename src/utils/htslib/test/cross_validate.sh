#!/bin/sh

#
# -----------------------------------------------------------------------------
# Author: James Bonfield.
#
# This cross validation script is designed to run the htslib test_view
# and cramtools.jar CRAM implementations to test compatibility between
# implementations.
#
# The test set may contain many dubious and ambiguous SAM cases, such as
# single base reads (is that quality "*" really meaning Q9 or no quality?).
# Some of these may fail one or the other implementations and be acceptable
# in the short-term, so to spot more important regressions we can tag
# specific cases as being known-pass or known-fail.
# -----------------------------------------------------------------------------
#

vers=3.0

cramtools_jar=$HOME/work/cram/cramtools/cramtools-$vers.jar

test_view="./test_view -o VERSION=$vers"

cramtools="/software/bin/java -Xmx4000m -jar $cramtools_jar"
cramtools="/software/bin/java -Xmx4000m -jar $cramtools_jar"

run_out() {
    out=$1; shift
    echo "$@ > $out"
    $@ > $out
}

run() {
    echo "$@"
    $@
}


sam_to_Ccram() {
    run_out _tmp.cram $test_view -C -t $1 $2
    #run_out _tmp.cram $HOME/io_lib/trunk/build.seq3/progs/scramble -r $1 -O CRAM $2
    if [ $? != 0 ]
    then
        crash=`expr $crash + 1`
        false
    fi
}

Ccram_to_sam() {
    run_out _tmp.sam $test_view -i REFERENCE=$1 _tmp.cram
    #run_out _tmp.sam $HOME/io_lib/trunk/build.seq3/progs/scramble -r $1 _tmp.cram

    if [ $? != 0 ]
    then
        crash=`expr $crash + 1`
        false
    fi
}

sam_to_Jcram() {
    run $cramtools cram -R $1 -I $2 -O _tmp.cram -n -Q --capture-all-tags
    if [ $? != 0 ]
    then
        crash=`expr $crash + 1`
        false
    fi
}

Jcram_to_sam() {
    run $cramtools bam -R $1 -I _tmp.cram -O _tmp.sam

    if [ $? != 0 ]
    then
        crash=`expr $crash + 1`
        false
    fi
}

compare_sam() {
    #run ./compare_sam.pl $i _tmp.sam -nomd -notemplate -unknownrg -Baux
    run ./compare_sam.pl $i _tmp.sam -nomd -Baux
    if [ $? != 0 ]
    then
        fails=`expr $fails + 1`
        false
    fi
}

trials=0
fails=0
crash=0

files=`ls -1 *#*.sam`

# Restrict to known workers from SAM->CRAM->CRAM in cramtools
#files="auxf#values.sam c1#bounds.sam c1#noseq.sam c1#pad1.sam c1#pad2.sam c1#pad3.sam c1#unknown.sam ce#1.sam ce#2.sam ce#5b.sam ce#large_seq.sam ce#tag_depadded.sam ce#tag_padded.sam ce#unmap.sam ce#unmap1.sam ce#unmap2.sam xx#large_aux.sam xx#large_aux2.sam xx#pair.sam xx#rg.sam xx#unsorted.sam"

for i in $files
do
    r=`echo $i | sed 's/#.*/.fa/'`
    echo "=== $i"

    # C to C
    trials=`expr $trials + 1`
    sam_to_Ccram $r $i && Ccram_to_sam $r && compare_sam $i _tmp.sam

    # Java to Java
    trials=`expr $trials + 1`
    sam_to_Jcram $r $i && Jcram_to_sam $r && compare_sam $i _tmp.sam

    # C to Java
    trials=`expr $trials + 1`
    sam_to_Ccram $r $i && Jcram_to_sam $r && compare_sam $i _tmp.sam

    # Java to C
    trials=`expr $trials + 1`
    sam_to_Jcram $r $i && Ccram_to_sam $r && compare_sam $i _tmp.sam
done

# Overcounts failures as an early fail can lead to 1 or 2 more fails.
echo ""
echo ============
echo No. tests: $trials
echo No. diffs: $fails
echo No. crash: $crash
