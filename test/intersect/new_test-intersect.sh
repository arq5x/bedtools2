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
#  Test intersection of a as bed from file vs b as bed from file
############################################################
echo "    intersect.new.t01...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a.bed -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bed from redirect vs b as bed from file
############################################################
echo "    intersect.new.t02...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bed from pipe vs b as bed from file
############################################################
echo "    intersect.new.t03...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a.bed | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bed from fifo vs b as bed from file
############################################################
echo "    intersect.new.t04...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a.bed) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as gzipped from file vs b as bed from file
############################################################
echo "    intersect.new.t05...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_gzipped.bed.gz -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as gzipped from redirect vs b as bed from file
############################################################
echo "    intersect.new.t06...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_gzipped.bed.gz > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as gzipped from pipe vs b as bed from file
############################################################
echo "    intersect.new.t07...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_gzipped.bed.gz | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as gzipped from fifo vs b as bed from file
############################################################
echo "    intersect.new.t08...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_gzipped.bed.gz) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bgzipped from file vs b as bed from file
############################################################
echo "    intersect.new.t09...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_bgzipped.bed.gz -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bgzipped from redirect vs b as bed from file
############################################################
echo "    intersect.new.t10...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_bgzipped.bed.gz > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bgzipped from pipe vs b as bed from file
############################################################
echo "    intersect.new.t11...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_bgzipped.bed.gz | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bgzipped from fifo vs b as bed from file
############################################################
echo "    intersect.new.t12...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_bgzipped.bed.gz) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bam from file vs b as bed from file
############################################################
echo "    intersect.new.t13...\c"
$BT intersect -a a.bam -b b.bed> obs
check obs aVSb.bam
rm obs


###########################################################
#  Test intersection of a as bam from redirect vs b as bed from file
############################################################
echo "    intersect.new.t14...\c"
$BT intersect -a - -b b.bed < a.bam> obs
check obs aVSb.bam
rm obs


###########################################################
#  Test intersection of a as bam from pipe vs b as bed from file
############################################################
echo "    intersect.new.t15...\c"
cat a.bam | $BT intersect -a - -b b.bed> obs
check obs aVSb.bam
rm obs


###########################################################
#  Test intersection of a as bam from fifo vs b as bed from file
############################################################
echo "    intersect.new.t16...\c"
$BT intersect -a <(cat a.bam) -b b.bed > obs
check obs aVSb.bam
rm obs


###########################################################
#  Test intersection of bam file containing both good reads
#  and those where both read and mate are unmapped vs b file
#  as bed.
############################################################
echo "    intersect.new.t17...\c"
echo \
"chr1	100	101	a2	255	-	100	200	0,0,0	1	100,	0,
chr1	100	110	a2	255	-	100	200	0,0,0	1	100,	0," > exp
$BT intersect -a a_with_bothUnmapped.bam -b b.bed -bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of bam file containing both good reads
#  and those where both read and mate are unmapped vs b file
#  as bed, with noHit (-v) option. 
############################################################
echo "    intersect.new.t18...\c"
echo \
"chr1	10	20	a1	255	+	10	20	0,0,0	1	10,	0,
.	-1	-1	FCC1MK2ACXX:1:1101:5780:51632#/1	0	.	-1	-1	-1	0,0,0	0	.	.
.	-1	-1	FCC1MK2ACXX:1:1101:5780:51632#/2	0	.	-1	-1	-1	0,0,0	0	.	.
.	-1	-1	FCC1MK2ACXX:1:1101:8137:99409#/1	0	.	-1	-1	-1	0,0,0	0	.	.
.	-1	-1	FCC1MK2ACXX:1:1101:8137:99409#/2	0	.	-1	-1	-1	0,0,0	0	.	.
.	-1	-1	FCC1MK2ACXX:1:1102:6799:2633#/1	0	.	-1	-1	-1	0,0,0	0	.	.
.	-1	-1	FCC1MK2ACXX:1:1102:6799:2633#/2	0	.	-1	-1	-1	0,0,0	0	.	." > exp
$BT intersect -a a_with_bothUnmapped.bam -b b.bed -bed -v > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of bam file containing read where query
#  is mapped and mate is not.
############################################################
echo "    intersect.new.t19...\c"
echo \
"chr1	98650	98704	FCC1MK2ACXX:1:1212:13841:9775#/1	0	+	98604	98704	0,0,0	1	100,	0," > exp
$BT intersect -a oneUnmapped.bam -b j1.bed -bed > obs
check obs exp
rm obs exp

###########################################################
#  Test intersection of bam file containing read where query
#  is mapped and mate is not, with noHit (-v) option.
############################################################
echo "    intersect.new.t20...\c"
echo \
"chr1	-1	-1	FCC1MK2ACXX:1:1212:13841:9775#/2	0	.	-1	-1	-1	0,0,0	0	.	." > exp
$BT intersect -a oneUnmapped.bam -b j1.bed -bed -v > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of bam file containing read where query
#  is unmapped but mate is mapped.
############################################################
echo "    intersect.new.t20.b...\c"
touch exp
$BT intersect -a queryUnmappedMateMappedCoordsInvalid.bam -b j1.bed -bed > obs
check obs exp
rm obs exp

###########################################################
#  Test intersection of bam file containing read where one
# mate is mapped and one is not, with noHit (-v) option.
############################################################
echo "    intersect.new.t20.c...\c"
echo ".	-1	-1	TTTACCTTT:4FSQ5P1:286:D2GA7ACXX:6:2316:20858:89646	0	.	-1	-1	-1	0,0,0	0	.	." > exp
$BT intersect -a queryUnmappedMateMappedCoordsInvalid.bam -b j1.bed -bed -v > obs
check obs exp
rm obs exp





###########################################################
#  Test intersection with -sorted, see that order specified
#  in genome file is enforced.
############################################################
echo "    intersect.new.t21...\c"
echo \
"Error: Sorted input specified, but the file chromOrderA.bed has the following record with a different sort order than the genomeFile human.hg19.genome
chr10	10	20	a3	100	+" > exp
$BT intersect -a chromOrderA.bed -b chromOrderB.bed -sorted -g human.hg19.genome 2>obs
check obs exp
rm obs exp


###########################################################
#  Test intersection with -sorted, see that chrom order
# change is ok as long as query and db have same order.
############################################################
echo "    intersect.new.t22...\c"
echo \
"chr1	15	20	a1	100	+
chr2	15	20	a2	100	+
chr10	15	20	a3	100	+
chr11	15	20	a4	100	+
chrX	15	20	a5	100	+
chrM	15	20	a6	100	+" > exp
$BT intersect -a chromOrderA.bed -b chromOrderB.bed -sorted > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection with -sorted, see that hits are correct
#  when sort order of files matches genome file sort order
############################################################
echo "    intersect.new.t23...\c"
echo \
"chr1	15	20	a1	100	+
chr2	15	20	a2	100	+
chr10	15	20	a3	100	+
chr11	15	20	a4	100	+
chrX	15	20	a5	100	+
chrM	15	20	a6	100	+" > exp
$BT intersect -a chromOrderA.bed -b chromOrderB.bed -sorted -g human.hg19.vSort.genome > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bed from file vs b as bed from file
############################################################
echo "    intersect.new.t24...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_withLargeHeader.bed -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bed from redirect vs b as bed from file
############################################################
echo "    intersect.new.t25...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_withLargeHeader.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bed from pipe vs b as bed from file
############################################################
echo "    intersect.new.t26...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_withLargeHeader.bed | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bed from fifo vs b as bed from file
############################################################
echo "    intersect.new.t27...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_withLargeHeader.bed) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as gzipped from file vs b as bed from file
############################################################
echo "    intersect.new.t28...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_withLargeHeader_gzipped.bed.gz -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as gzipped from redirect vs b as bed from file
############################################################
echo "    intersect.new.t29...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_withLargeHeader_gzipped.bed.gz > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as gzipped from pipe vs b as bed from file
############################################################
echo "    intersect.new.t30...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_withLargeHeader_gzipped.bed.gz | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as gzipped from fifo vs b as bed from file
############################################################
echo "    intersect.new.t31...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_withLargeHeader_gzipped.bed.gz) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bgzipped from file vs b as bed from file
############################################################
echo "    intersect.new.t32...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_withLargeHeader_bgzipped.bed.gz -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bgzipped from redirect vs b as bed from file
############################################################
echo "    intersect.new.t33...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_withLargeHeader_bgzipped.bed.gz > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bgzipped from pipe vs b as bed from file
############################################################
echo "    intersect.new.t34...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_withLargeHeader_bgzipped.bed.gz | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bgzipped from fifo vs b as bed from file
############################################################
echo "    intersect.new.t35...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_withLargeHeader_bgzipped.bed.gz) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bed from file 
#  vs b as bed from file, and print header
############################################################
echo "    intersect.new.t36...\c"
$BT intersect -a a_withLargeHeader.bed -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bed from 
#  redirect vs b as bed from file, and print header
############################################################
echo "    intersect.new.t37...\c"
$BT intersect -a - -b b.bed < a_withLargeHeader.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bed from pipe
#  vs b as bed from file, and print header
############################################################
echo "    intersect.new.t38...\c"
cat a_withLargeHeader.bed | $BT intersect -a - -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bed from fifo
#  vs b as bed from file, and print header
############################################################
echo "    intersect.new.t39...\c"
$BT intersect -a <(cat a_withLargeHeader.bed) -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as gzipped from
#  file vs b as bed from file, and print header
############################################################
echo "    intersect.new.t40...\c"
$BT intersect -a a_withLargeHeader_gzipped.bed.gz -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as gzipped from
#  redirect vs b as bed from file, and print header
############################################################
echo "    intersect.new.t41...\c"
$BT intersect -a - -b b.bed < a_withLargeHeader_gzipped.bed.gz -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as gzipped from
#  pipe vs b as bed from file, and print header
############################################################
echo "    intersect.new.t42...\c"
cat a_withLargeHeader_gzipped.bed.gz | $BT intersect -a - -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as gzipped from
#  fifo vs b as bed from file, and print header
############################################################
echo "    intersect.new.t43...\c"
$BT intersect -a <(cat a_withLargeHeader_gzipped.bed.gz) -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bgzipped from
#  file vs b as bed from file, and print header
############################################################
echo "    intersect.new.t44...\c"
$BT intersect -a a_withLargeHeader_bgzipped.bed.gz -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bgzipped from
#  redirect vs b as bed from file, and print header
############################################################
echo "    intersect.new.t45...\c"
$BT intersect -a - -b b.bed < a_withLargeHeader_bgzipped.bed.gz -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bgzipped from
#  pipe vs b as bed from file, and print header
############################################################
echo "    intersect.new.t46...\c"
cat a_withLargeHeader_bgzipped.bed.gz | $BT intersect -a - -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bgzipped from
#  fifo vs b as bed from file, and print header
############################################################
echo "    intersect.new.t47...\c"
$BT intersect -a <(cat a_withLargeHeader_bgzipped.bed.gz) -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test enforcement of sort order when using -sorted option.
#  Show that detect chrom "jumping", i.e. chr1, chr2, then
#  back to chr1
############################################################
echo "    intersect.new.t48...\c"
echo "Error: Sorted input specified" > exp
$BT intersect -a chromsOutOfOrder.bed -b b.bed -sorted 2>&1 > /dev/null | grep "Error: Sorted input specified, " | cut -f1 -d ',' > obs
check obs exp
rm obs exp

###########################################################
#  Test enforcement of sort order when using -sorted option.
#  Show that detect chrom "jumping", i.e. chr1, chr2, then
#  back to chr1
############################################################
echo "    intersect.new.t49...\c"
echo "Error: Sorted input specified" > exp
$BT intersect -a recordsOutOfOrder.bed -b b.bed -sorted 2>&1 > /dev/null | grep "Error: Sorted input specified, " | cut -f1 -d ',' > obs
check obs exp
rm obs exp

###########################################################
#  Test that we throw an error for non-existant files
############################################################
echo "    intersect.new.t50...\c"
echo "Error: Unable to open file nonExistantFile.bed. Exiting." > exp
$BT intersect -a nonExistantFile.bed -b b.bed -sorted 2>&1 > /dev/null | cat - >obs
check obs exp
rm obs exp

###########################################################
#  Test that we allow empty query files
############################################################
echo "    intersect.new.t51...\c"
touch exp
touch dummy.txt
$BT intersect -a dummy.txt -b b.bed 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs exp dummy.txt


###########################################################
#  Test that we throw an error for unrecognized arguments
############################################################
echo "    intersect.new.t52...\c"
echo "***** ERROR: Unrecognized parameter: -wrongArg *****" > exp
$BT intersect -a a.bed -b b.bed -wrongArg 2>&1 > /dev/null | head -2 | tail -1 > obs
check obs exp
rm obs exp


###########################################################
#  Test that we can process a Bam file with no text in 
#  it's header. 
############################################################
echo "    intersect.new.t53...\c"
$BT intersect -a gdc.bam -b gdc.bam -bed > obs
check obs gdc_exp
rm obs

###########################################################
#  Test that if the query is BAM, and bed output is not 
#  explicit, and they asked for an output option that is not
# valid with BAM output, an error is thrown.
############################################################
echo "    intersect.new.t54...\c"
echo "***** ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option. *****" > exp
$BT intersect -a a.bam -b b.bed  -wao 2>&1 > /dev/null  | head -2 | tail -1 > obs
check obs exp
rm obs exp


###########################################################
# Test that gff files work correctly
############################################################
echo "    intersect.new.t55...\c"
echo \
"chr2L	.	UTR	51	70	0	+	.	ID=mRNA:xs2:UTR:41-70;Parent=mRNA:xs2;
chr2L	.	CDS	71	100	0	+	.	ID=mRNA:xs2:CDS:71-130;Parent=mRNA:xs2;
chr2L	.	exon	51	100	0	+	.	ID=mRNA:xs2:exon:41-130;Parent=mRNA:xs2;
chr2L	.	mRNA	51	100	0	+	.	ID=mRNA:xs2;Parent=g2;
chr2L	.	gene	51	100	0	+	.	ID=g2;" > exp
$BT intersect -a gdc.gff -b gdc_one.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test that we allow empty database files, unsorted
############################################################
echo "    intersect.new.t56...\c"
touch exp
touch dummy.txt
$BT intersect -a a.bed -b dummy.txt 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs exp dummy.txt

###########################################################
#  Test that we allow empty database files, sorted
############################################################
echo "    intersect.new.t57...\c"
touch exp
touch dummy.txt
$BT intersect -a a.bed -b dummy.txt -sorted 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs exp dummy.txt


###########################################################
#  Test that we allow empty database files, unsorted, with
# -v (noHit) option
############################################################
echo "    intersect.new.t58...\c"
touch dummy.txt
$BT intersect -a a.bed -b dummy.txt -v > obs
check obs a.bed
rm obs dummy.txt


###########################################################
#  Test that an empty query with header run with -header
# option will print header
############################################################
echo "    intersect.new.t59...\c"
echo "#Random Header" >dummy.txt
$BT intersect -a dummy.txt -b a.bed -header > obs
check obs dummy.txt
rm obs dummy.txt


###########################################################
#  Test that an empty query with header, gzipped, that
#  runs with -header option will print header
############################################################
echo "    intersect.new.t60...\c"
echo "#Random Header" >dummy.txt
gzip dummy.txt
echo "#Random Header" >exp
$BT intersect -a dummy.txt.gz -b a.bed -header > obs
check obs exp
rm obs dummy.txt.gz exp

###########################################################
#  Test that an empty query with header, bgzipped, that
#  runs with -header option will print header
############################################################
echo "    intersect.new.t61...\c"
echo "#Random Header" >dummy.txt
bgzip dummy.txt
echo "#Random Header" >exp
$BT intersect -a dummy.txt.gz -b a.bed -header > obs
check obs exp
rm obs dummy.txt.gz exp


###########################################################
#  Test that an empty VCF query with header that
#  runs with -header option will print header
############################################################
echo "    intersect.new.t62...\c"
$BT intersect -a headerOnly.vcf -b a.bed -header > obs
check obs headerOnly.vcf
rm obs


###########################################################
#  Test that files with DOS newline characters, '\r',
#  and/or extra tabs at end of line are handled
############################################################
echo "    intersect.new.t63...\c"
echo \
"Error: Type checker found wrong number of fields while tokenizing data line."   >exp
$BT intersect -a dosLineChar_a.bed -b dosLineCharWithExtraTab_b.bed -v 2>&1 > /dev/null | cat - > obs
check exp obs
rm exp obs

###########################################################
#  Test that files with no newlines at all are handled
############################################################
echo "    intersect.new.t64...\c"
echo "chr17	7577068	7577157" > exp
$BT intersect -a oneRecordNoNewline.bed -b oneRecordNoNewline.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test zero length intersections in non-bam files.
############################################################
echo "    intersect.new.t65...\c"
echo \
"chr1	5	15	r1	chr1	9	9	m3	0
chr1	7	12	r3	chr1	9	9	m3	0" > exp
$BT intersect -a a_testZeroLen.bed -b b_testZeroLen.bed -wo > obs
check exp obs
rm exp obs

###########################################################
#  Test zero length intersections in non-bam files, -sorted
############################################################
echo "    intersect.new.t66...\c"
echo \
"chr1	5	15	r1	chr1	9	9	m3	0
chr1	7	12	r3	chr1	9	9	m3	0" > exp
$BT intersect -a a_testZeroLen.bed -b b_testZeroLen.bed -wo -sorted> obs
check exp obs
rm exp obs

###########################################################
#  Test vcf struct var intersection
############################################################
echo "    intersect.new.t67a...\c"
echo \
"19	252806	791255	G	<DEL>	70.90	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-389,-4611;END=253195;STR=+-:4;IMPRECISE;CIPOS=-2,137;CIEND=0,0;EVENT=791255;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	256900	791255	G	T	70.90	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-389,-4611;END=253195;STR=+-:4;IMPRECISE;CIPOS=-2,137;CIEND=0,0;EVENT=791255;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	260365	791256	C	<DEL>	33.71	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-680;END=261045;STR=+-:4;IMPRECISE;CIPOS=-1,257;CIEND=0,0;EVENT=791256;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=upstream_gene_variant|||ENSG00000271846|CTD-3113P16.9|ENST00000607399|||||processed_pseudogene	19	260800	791256	C	<INS>	33.71	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=680;END=261045;STR=+-:4;IMPRECISE;CIPOS=-1,257;CIEND=0,0;EVENT=791256;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=upstream_gene_variant|||ENSG00000271846|CTD-3113P16.9|ENST00000607399|||||processed_pseudogene
19	265134	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	265500	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	265986	791258	A	<DEL>	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	265500	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	265986	791258	A	<DEL>	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	266003	791258	A	C	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||" > exp
$BT intersect -a a_vcfSVtest.vcf -b b_vcfSVtest.vcf -wa -wb >obs
check exp obs
rm exp obs

echo "    intersect.new.t67b...\c"
echo -n "" > exp
$BT intersect -a a_vcfSVtest.vcf -b b_vcfSVtest.vcf -wa -wb -v >obs
check exp obs
rm exp obs

###########################################################
#  Test vcf struct var intersection, sorted
############################################################
echo "    intersect.new.t68a...\c"
echo \
"19	252806	791255	G	<DEL>	70.90	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-389,-4611;END=253195;STR=+-:4;IMPRECISE;CIPOS=-2,137;CIEND=0,0;EVENT=791255;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	256900	791255	G	T	70.90	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-389,-4611;END=253195;STR=+-:4;IMPRECISE;CIPOS=-2,137;CIEND=0,0;EVENT=791255;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	260365	791256	C	<DEL>	33.71	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-680;END=261045;STR=+-:4;IMPRECISE;CIPOS=-1,257;CIEND=0,0;EVENT=791256;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=upstream_gene_variant|||ENSG00000271846|CTD-3113P16.9|ENST00000607399|||||processed_pseudogene	19	260800	791256	C	<INS>	33.71	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=680;END=261045;STR=+-:4;IMPRECISE;CIPOS=-1,257;CIEND=0,0;EVENT=791256;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=upstream_gene_variant|||ENSG00000271846|CTD-3113P16.9|ENST00000607399|||||processed_pseudogene
19	265134	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	265500	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	265986	791258	A	<DEL>	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	265500	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	265986	791258	A	<DEL>	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	266003	791258	A	C	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||" > exp
$BT intersect -a a_vcfSVtest.vcf -b b_vcfSVtest.vcf -wa -wb -sorted >obs
check exp obs
rm exp obs

echo "    intersect.new.t68b...\c"
echo -n "" > exp
$BT intersect -a a_vcfSVtest.vcf -b b_vcfSVtest.vcf -wa -wb -v -sorted>obs
check exp obs
rm exp obs

###########################################################
#  Test intersect -loj with multiple databases
############################################################
echo "    intersect.new.t69...\c"
echo \
"1	100	200	a1	ax	1	1	100	200	b1	bx
1	300	400	a2	ay	.	.	-1	-1	.	.
1	400	500	a3	az	.	.	-1	-1	.	.
2	500	600	a4	aq	.	.	-1	-1	.	." > exp
$BT intersect -a null_a.bed -b null_b.bed null_c.bed -loj > obs
check exp obs
rm exp obs

###########################################################
#  Test intersect -loj with multiple databases and -names
############################################################
echo "    intersect.new.t70...\c"
echo \
"1	100	200	a1	ax	b	1	100	200	b1	bx
1	300	400	a2	ay	.	.	-1	-1	.	.
1	400	500	a3	az	.	.	-1	-1	.	.
2	500	600	a4	aq	.	.	-1	-1	.	." > exp
$BT intersect -a null_a.bed -b null_b.bed null_c.bed -loj -names b c > obs
check exp obs
rm exp obs

###########################################################
#  Test intersect -wao with multiple databases
############################################################
echo "    intersect.new.t71...\c"
echo \
"1	100	200	a1	ax	1	1	100	200	b1	bx	100
1	300	400	a2	ay	.	.	-1	-1	.	.	0
1	400	500	a3	az	.	.	-1	-1	.	.	0
2	500	600	a4	aq	.	.	-1	-1	.	.	0" > exp
$BT intersect -a null_a.bed -b null_b.bed null_c.bed -wao > obs
check exp obs
rm exp obs

###########################################################
#  Test intersect -wao with multiple databases and -names
############################################################
echo "    intersect.new.t72...\c"
echo \
"1	100	200	a1	ax	b	1	100	200	b1	bx	100
1	300	400	a2	ay	.	.	-1	-1	.	.	0
1	400	500	a3	az	.	.	-1	-1	.	.	0
2	500	600	a4	aq	.	.	-1	-1	.	.	0" > exp
$BT intersect -a null_a.bed -b null_b.bed null_c.bed -wao -names b c > obs
check exp obs
rm exp obs
