BT=${BT-../../bin/bedtools}

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
samtools view -Sb pair-sorted.sam > pair-sorted.bam 2>/dev/null
samtools view -Sb pair-unsorted.sam > pair-unsorted.bam 2>/dev/null
samtools view -Sb pair-and-missingmate.sam > pair-and-missingmate.bam 2>/dev/null
samtools view -Sb missingmate-and-pair.sam > missingmate-and-pair.bam 2>/dev/null
samtools view -Sb pair-with-name-containing-slash.sam > pair-with-name-containing-slash.bam 2>/dev/null
samtools view -Sb proper_pair_plus_strand.sam > proper_pair_plus_strand.bam 2> /dev/null
samtools view -Sb proper_pair_plus_strand.sam > proper_pair_plus_strand.bam 2> /dev/null
samtools view -Sb proper_pair_minus_strand.sam > proper_pair_minus_strand.bam 2> /dev/null
samtools view -Sb proper_pair_minus_strand.sam > proper_pair_minus_strand.bam 2> /dev/null
samtools view -Sb proper_pair_spliced.sam > proper_pair_spliced.bam 2> /dev/null
samtools view -Sb proper_pair_spliced.sam > proper_pair_spliced.bam 2> /dev/null
samtools view -Sb proper_pair_overlap.sam > proper_pair_overlap.bam 2> /dev/null
samtools view -Sb proper_pair_overlap.sam > proper_pair_overlap.bam 2> /dev/null
samtools view -Sb proper_pair_bad_mapq.sam > proper_pair_bad_mapq.bam 2> /dev/null
samtools view -Sb proper_pair_bad_mapq.sam > proper_pair_bad_mapq.bam 2> /dev/null
samtools view -Sb not_proper_pair.sam > not_proper_pair.bam 2> /dev/null
samtools view -Sb not_proper_pair.sam > not_proper_pair.bam 2> /dev/null
samtools view -Sb bug_proper_pair_different_chrom.sam > bug_proper_pair_different_chrom.bam 2> /dev/null
samtools view -Sb bug_proper_pair_different_chrom.sam > bug_proper_pair_different_chrom.bam 2> /dev/null



##################################################################
#  Test paired reads when sorted by name and position
##################################################################
echo "    pairedbamtobed12.t1...\c"
echo \
"chr14	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162" > exp
$BT pairedbamtobed12 -i pair-sorted.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test paired reads when sorted by name and not by position
##################################################################
echo "    pairedbamtobed12.t2...\c"
echo \
"chr14	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162" > exp
$BT pairedbamtobed12 -i pair-unsorted.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test skipping missing mate first read
##################################################################
echo "    pairedbamtobed12.t3...\c"
echo \
"chr4	8210431	8210761	M00528:19:000000000-A88YD:1:1101:2318:12845	120	+	8210431	8210458	255,0,0	2	27,21	0,309" > exp

$BT pairedbamtobed12 -i missingmate-and-pair.bam > obs 2> /dev/null
check obs exp
rm obs exp


##################################################################
#  Test skipping missing mate first read (test stderr)
##################################################################
echo "    pairedbamtobed12.t3.stderr...\c"
echo \
"*****WARNING: Query M00528:19:000000000-A88YD:1:1101:2241:12366 is not followed by his mate in your BAM file. Skipping" > exp
$BT pairedbamtobed12 -i missingmate-and-pair.bam > /dev/null 2> obs
check obs exp
rm obs exp


##################################################################
#  Test skipping missing mate last read
##################################################################
echo "    pairedbamtobed12.t4...\c"
echo \
"chr14	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162" > exp
$BT pairedbamtobed12 -i pair-and-missingmate.bam > obs 2> /dev/null
check obs exp
rm obs exp


##################################################################
#  Test skipping missing mate last read (stderr)
##################################################################
echo "    pairedbamtobed12.t4.stderr...\c"
echo \
"*****WARNING: Query M00528:19:000000000-A88YD:1:1101:2318:12845 is the last read and has no mate. Skip and exit. " > exp
$BT pairedbamtobed12 -i pair-and-missingmate.bam > /dev/null 2> obs
check obs exp
rm obs exp
#exit

##################################################################
#  Test paired reads with name ending with '/*'
##################################################################
echo "    pairedbamtobed12.t5...\c"
echo \
"chr14	50053297	50053480	M00528:19:000000000-A88YD:1:1101:2241:12366	0	+	50053297	50053324	255,0,0	2	27,21	0,162" > exp
$BT pairedbamtobed12 -i pair-with-name-containing-slash.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test proper pair on the plus strand
##################################################################
echo "    pairedbamtobed12.t6...\c"
echo \
"chr1	0	99	proper_pair_plus_strand	80	+	0	30	255,0,0	2	30,30	0,69" > exp
$BT pairedbamtobed12 -i proper_pair_plus_strand.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test proper pair on the minus strand
##################################################################
echo "    pairedbamtobed12.t7...\c"
echo \
"chr1	0	99	proper_pair_minus_strand	80	-	69	99	255,0,0	2	30,30	0,69" > exp
$BT pairedbamtobed12 -i proper_pair_minus_strand.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test proper pair spliced
##################################################################
echo "    pairedbamtobed12.t8...\c"
echo \
"chr1	0	99	proper_pair_spliced	80	+	0	40	255,0,0	3	20,10,30	0,30,69" > exp
$BT pairedbamtobed12 -i proper_pair_spliced.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test proper pair overlap
##################################################################
echo "    pairedbamtobed12.t9...\c"
echo \
"chr1	0	49	proper_pair_overlap	80	+	0	30	255,0,0	2	30,30	0,19" > exp
$BT pairedbamtobed12 -i proper_pair_overlap.bam > obs
check obs exp
rm obs exp


##################################################################
#  Test skipping proper pair with bad mapq when using -qual argument
##################################################################
echo "    pairedbamtobed12.t10...\c"
touch exp
$BT pairedbamtobed12 -i proper_pair_bad_mapq.bam -qual 11 > obs
check obs exp
rm obs exp


##################################################################
#  Test skip not proper pair
##################################################################
echo "    pairedbamtobed12.t11...\c"
touch exp
$BT pairedbamtobed12 -i not_proper_pair.bam > obs 2> /dev/null
check obs exp
rm obs exp

##################################################################
#  Test skip not proper pair (stderr)
##################################################################
echo "    pairedbamtobed12.t11.stderr...\c"
echo \
"*****WARNING: Query not_proper_pair is not followed by his mate in your BAM file. Skipping
*****WARNING: Query not_proper_pair is the last read and has no mate. Skip and exit. " > exp
$BT pairedbamtobed12 -i not_proper_pair.bam > /dev/null 2> obs 
check obs exp
rm obs exp

##################################################################
#  Test skip proper when reads are on different chromosomes
##################################################################
echo "    pairedbamtobed12.t12...\c"
touch exp
$BT pairedbamtobed12 -i bug_proper_pair_different_chrom.bam > obs 2> /dev/null
check obs exp
rm obs exp

##################################################################
#  Test skip proper when reads are on different chromosomes (stderr)
##################################################################
echo "    pairedbamtobed12.t12.stderr...\c"
echo \
"*****WARNING: Query bug_proper_pair_different_chrom is not on the same chromosome than his mate. Skipping
*****WARNING: Query bug_proper_pair_different_chrom is the last read and has no mate. Skip and exit. " > exp
$BT pairedbamtobed12 -i bug_proper_pair_different_chrom.bam > /dev/null 2> obs
check obs exp
rm obs exp


rm *.bam
