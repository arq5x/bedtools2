#!/bin/sh
# The training run for profile-guided optimization (PGO in the Makefile).
#
#   sh scripts/pgo-train.sh <instrumented bedtools> <scratch directory>
#
# Run from the top of the source tree. The profile weights the code by how
# often these calls run it: the large files decide which paths are hot, and the
# small ones keep the other modes and formats from being laid out as cold.
set -e
BT=$1
S=$2
T=test/intersect/sortAndNaming/bigTests
I=test/intersect
A=$T/q500K.bed
D=$T/db500K.bed
mkdir -p $S

# 500,000 BED6 records against 500,000 others, about 3% of A overlapping.
for m in "" "-wa -wb" "-u" "-c" "-v" "-sorted"; do
	$BT intersect -a $A -b $D $m > /dev/null
done
$BT subtract -a $A -b $D > /dev/null
$BT coverage -a $A -b $D > /dev/null

# Every A record overlapping: A as BED6+ (column 6 not a strand) against
# copies of its own intervals as BED6+4 and as BED3.
awk 'BEGIN { OFS = "\t" } { $6 = "ELS"; print }' $A > $S/aplus.bed
awk 'BEGIN { OFS = "\t" } { print $0, "x", 7, 8, "y" }' $A > $S/b64.bed
cut -f1-3 $A > $S/b3.bed
for m in "" "-wa -wb" "-u"; do
	$BT intersect -a $S/aplus.bed -b $S/b64.bed $m > /dev/null
done
$BT intersect -a $S/aplus.bed -b $S/b3.bed -c > /dev/null

# BED12 with two blocks per record under -split, gzip-compressed when gzip is
# available.
awk 'BEGIN { OFS = "\t" } { print $1, $2, $3, $4, $5, $6, $2, $3, 0, 2, "40,40,", "0,60," }' $A > $S/a12.bed
if gzip -c $S/a12.bed > $S/a12.bed.gz && gzip -c $S/b3.bed > $S/b3.bed.gz; then
	$BT intersect -a $S/a12.bed.gz -b $S/b3.bed.gz -split -c > /dev/null
else
	$BT intersect -a $S/a12.bed -b $S/b3.bed -split -c > /dev/null
fi

# Small inputs for the other modes and formats.
for m in "" -u -v -c -C "-wa -wb" -wo -wao -loj -s -S "-f 0.5 -r" "-f 0.5 -F 0.5 -e" -sortout -sorted -split; do
	$BT intersect -a $I/a.bed -b $I/b.bed $m > /dev/null
done
$BT intersect -a $I/a.bed -b $I/b.bed $I/c.bed -names b c -wa -wb > /dev/null
$BT intersect -a $I/a.bam -b $I/b.bed > /dev/null
$BT intersect -a $I/a.bam -b $I/b.bed -bed > /dev/null
$BT intersect -a $I/a_bgzipped.bed.gz -b $I/b.bed > /dev/null
$BT intersect -a $I/a_gzipped.bed.gz -b $I/b.bed -c > /dev/null
$BT intersect -a $I/a_withLargeHeader_bgzipped.bed.gz -b $I/b.bed -header > /dev/null
$BT intersect -a $I/a_vcfSVtest.vcf -b $I/b_vcfSVtest.vcf > /dev/null
$BT intersect -a $I/gdc.gff -b $I/gdc_one.bed > /dev/null
$BT intersect -a $I/blocks.bed12 -b $I/three_blocks_match.bed -split > /dev/null

# Other common subcommands, each on files of the same size.
$BT sort -i $D > /dev/null
$BT merge -i $D -c 5 -o count > /dev/null
$BT closest -a $A -b $S/b3.bed -d -io > /dev/null
$BT map -a $A -b $S/b64.bed -c 5 -o mean > /dev/null
$BT window -a $S/aplus.bed -b $S/b3.bed -w 500 -c > /dev/null
$BT jaccard -a $S/aplus.bed -b $S/b3.bed > /dev/null

# Others, on the chr1 records of the same files.
Q1=$T/qchr1.bed
D1=$T/dbchr1.bed
G=$I/human.hg19.genome
$BT multiinter -i $Q1 $D1 > /dev/null
$BT cluster -i $Q1 > /dev/null
$BT groupby -i $Q1 -g 1 -c 5 -o sum > /dev/null
$BT spacing -i $Q1 > /dev/null
$BT slop -i $Q1 -g $G -b 100 > /dev/null
$BT flank -i $Q1 -g $G -b 100 > /dev/null
$BT complement -i $Q1 -g $G > /dev/null
$BT shuffle -i $Q1 -g $G -seed 1 > /dev/null
