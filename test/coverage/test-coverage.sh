BT=../../bin/bedtools

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}

###########################################################
###########################################################
#                       BAM files                         #
###########################################################
###########################################################
samtools view -Sb one_block.sam > one_block.bam 2>/dev/null
samtools view -Sb two_blocks.sam > two_blocks.bam 2>/dev/null
samtools view -Sb three_blocks.sam > three_blocks.bam 2>/dev/null
samtools view -Sb sam-w-del.sam > sam-w-del.bam 2>/dev/null



##################################################################
#  Test three blocks without -split
##################################################################
# echo "    coverage.t1...\c"
# echo \
# "chr1	0	50	1" > exp
# $BT coverage -abam three_blocks.bam -b three_blocks_nomatch.bed > obs
# check obs exp
# rm obs exp


rm *.bam