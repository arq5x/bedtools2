BT=${BT-../../bin/bedtools}
htsutil=${htsutil-../htsutil}

FAILURES=0;

check()
{
	if diff $1 $2; then
    	echo ok
	else
    	FAILURES=$(expr $FAILURES + 1);
		echo fail
	fi
}

bam_check() 
{
	if diff <($htsutil viewbamrecords $1) <($htsutil viewbamrecords $2)
	then
		echo ok
	else
		FAILURES$(expr $FAILURES + 1);
		echo fail
	fi
}

###########################################################
#  Test intersection of a as bed from file vs b as bed from file
############################################################
echo -e "    intersect.new.t01...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a.bed -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bed from redirect vs b as bed from file
############################################################
echo -e "    intersect.new.t02...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bed from pipe vs b as bed from file
############################################################
echo -e "    intersect.new.t03...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a.bed | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bed from fifo vs b as bed from file
############################################################
echo -e "    intersect.new.t04...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a.bed) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as gzipped from file vs b as bed from file
############################################################
echo -e "    intersect.new.t05...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_gzipped.bed.gz -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as gzipped from redirect vs b as bed from file
############################################################
echo -e "    intersect.new.t06...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_gzipped.bed.gz > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as gzipped from pipe vs b as bed from file
############################################################
echo -e "    intersect.new.t07...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_gzipped.bed.gz | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as gzipped from fifo vs b as bed from file
############################################################
echo -e "    intersect.new.t08...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_gzipped.bed.gz) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bgzipped from file vs b as bed from file
############################################################
echo -e "    intersect.new.t09...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_bgzipped.bed.gz -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bgzipped from redirect vs b as bed from file
############################################################
echo -e "    intersect.new.t10...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_bgzipped.bed.gz > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bgzipped from pipe vs b as bed from file
############################################################
echo -e "    intersect.new.t11...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_bgzipped.bed.gz | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bgzipped from fifo vs b as bed from file
############################################################
echo -e "    intersect.new.t12...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_bgzipped.bed.gz) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a as bam from file vs b as bed from file
############################################################
echo -e "    intersect.new.t13...\c"
$BT intersect -a a.bam -b b.bed> obs
bam_check obs aVSb.bam
rm obs


###########################################################
#  Test intersection of a as bam from redirect vs b as bed from file
############################################################
echo -e "    intersect.new.t14...\c"
$BT intersect -a - -b b.bed < a.bam> obs
bam_check obs aVSb.bam
rm obs


###########################################################
#  Test intersection of a as bam from pipe vs b as bed from file
############################################################
echo -e "    intersect.new.t15...\c"
cat a.bam | $BT intersect -a - -b b.bed> obs
bam_check obs aVSb.bam
rm obs


###########################################################
#  Test intersection of a as bam from fifo vs b as bed from file
############################################################
echo -e "    intersect.new.t16...\c"
$BT intersect -a <(cat a.bam) -b b.bed > obs
bam_check obs aVSb.bam
rm obs


###########################################################
#  Test intersection of bam file containing both good reads
#  and those where both read and mate are unmapped vs b file
#  as bed.
############################################################
echo -e "    intersect.new.t17...\c"
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
echo -e "    intersect.new.t18...\c"
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
echo -e "    intersect.new.t19...\c"
echo \
"chr1	98650	98704	FCC1MK2ACXX:1:1212:13841:9775#/1	0	+	98604	98704	0,0,0	1	100,	0," > exp
$BT intersect -a oneUnmapped.bam -b j1.bed -bed > obs
check obs exp
rm obs exp

###########################################################
#  Test intersection of bam file containing read where query
#  is mapped and mate is not, with noHit (-v) option.
############################################################
echo -e "    intersect.new.t20...\c"
echo \
"chr1	-1	-1	FCC1MK2ACXX:1:1212:13841:9775#/2	0	.	-1	-1	-1	0,0,0	0	.	." > exp
$BT intersect -a oneUnmapped.bam -b j1.bed -bed -v > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of bam file containing read where query
#  is unmapped but mate is mapped.
############################################################
echo -e "    intersect.new.t20.b...\c"
touch exp
$BT intersect -a queryUnmappedMateMappedCoordsInvalid.bam -b j1.bed -bed > obs
check obs exp
rm obs exp

###########################################################
#  Test intersection of bam file containing read where one
# mate is mapped and one is not, with noHit (-v) option.
############################################################
echo -e "    intersect.new.t20.c...\c"
echo ".	-1	-1	TTTACCTTT:4FSQ5P1:286:D2GA7ACXX:6:2316:20858:89646	0	.	-1	-1	-1	0,0,0	0	.	." > exp
$BT intersect -a queryUnmappedMateMappedCoordsInvalid.bam -b j1.bed -bed -v > obs
check obs exp
rm obs exp





###########################################################
#  Test intersection with -sorted, see that order specified
#  in genome file is enforced.
############################################################
echo -e "    intersect.new.t21...\c"
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
echo -e "    intersect.new.t22...\c"
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
echo -e "    intersect.new.t23...\c"
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
echo -e "    intersect.new.t24...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_withLargeHeader.bed -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bed from redirect vs b as bed from file
############################################################
echo -e "    intersect.new.t25...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_withLargeHeader.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bed from pipe vs b as bed from file
############################################################
echo -e "    intersect.new.t26...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_withLargeHeader.bed | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bed from fifo vs b as bed from file
############################################################
echo -e "    intersect.new.t27...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_withLargeHeader.bed) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as gzipped from file vs b as bed from file
############################################################
echo -e "    intersect.new.t28...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_withLargeHeader_gzipped.bed.gz -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as gzipped from redirect vs b as bed from file
############################################################
echo -e "    intersect.new.t29...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_withLargeHeader_gzipped.bed.gz > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as gzipped from pipe vs b as bed from file
############################################################
echo -e "    intersect.new.t30...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_withLargeHeader_gzipped.bed.gz | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as gzipped from fifo vs b as bed from file
############################################################
echo -e "    intersect.new.t31...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a <(cat a_withLargeHeader_gzipped.bed.gz) -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bgzipped from file vs b as bed from file
############################################################
echo -e "    intersect.new.t32...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a a_withLargeHeader_bgzipped.bed.gz -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bgzipped from redirect vs b as bed from file
############################################################
echo -e "    intersect.new.t33...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
$BT intersect -a - -b b.bed < a_withLargeHeader_bgzipped.bed.gz > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bgzipped from pipe vs b as bed from file
############################################################
echo -e "    intersect.new.t34...\c"
echo \
"chr1	100	101	a2	2	-
chr1	100	110	a2	2	-" > exp
cat a_withLargeHeader_bgzipped.bed.gz | $BT intersect -a - -b b.bed > obs
check obs exp
rm obs exp


###########################################################
#  Test intersection of a with large header as bgzipped from fifo vs b as bed from file
############################################################
echo -e "    intersect.new.t35...\c"
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
echo -e "    intersect.new.t36...\c"
$BT intersect -a a_withLargeHeader.bed -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bed from 
#  redirect vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t37...\c"
$BT intersect -a - -b b.bed < a_withLargeHeader.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bed from pipe
#  vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t38...\c"
cat a_withLargeHeader.bed | $BT intersect -a - -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bed from fifo
#  vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t39...\c"
$BT intersect -a <(cat a_withLargeHeader.bed) -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as gzipped from
#  file vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t40...\c"
$BT intersect -a a_withLargeHeader_gzipped.bed.gz -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as gzipped from
#  redirect vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t41...\c"
$BT intersect -a - -b b.bed < a_withLargeHeader_gzipped.bed.gz -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as gzipped from
#  pipe vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t42...\c"
cat a_withLargeHeader_gzipped.bed.gz | $BT intersect -a - -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as gzipped from
#  fifo vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t43...\c"
$BT intersect -a <(cat a_withLargeHeader_gzipped.bed.gz) -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bgzipped from
#  file vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t44...\c"
$BT intersect -a a_withLargeHeader_bgzipped.bed.gz -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bgzipped from
#  redirect vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t45...\c"
$BT intersect -a - -b b.bed < a_withLargeHeader_bgzipped.bed.gz -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bgzipped from
#  pipe vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t46...\c"
cat a_withLargeHeader_bgzipped.bed.gz | $BT intersect -a - -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test intersection of a with large header as bgzipped from
#  fifo vs b as bed from file, and print header
############################################################
echo -e "    intersect.new.t47...\c"
$BT intersect -a <(cat a_withLargeHeader_bgzipped.bed.gz) -b b.bed -header > obs
check obs aWithHeaderVsB.txt
rm obs


###########################################################
#  Test enforcement of sort order when using -sorted option.
#  Show that detect chrom "jumping", i.e. chr1, chr2, then
#  back to chr1
############################################################
echo -e "    intersect.new.t48...\c"
echo "Error: Sorted input specified" > exp
$BT intersect -a chromsOutOfOrder.bed -b b.bed -sorted 2>&1 > /dev/null | grep "Error: Sorted input specified, " | cut -f1 -d ',' > obs
check obs exp
rm obs exp

###########################################################
#  Test enforcement of sort order when using -sorted option.
#  Show that detect chrom "jumping", i.e. chr1, chr2, then
#  back to chr1
############################################################
echo -e "    intersect.new.t49...\c"
echo "Error: Sorted input specified" > exp
$BT intersect -a recordsOutOfOrder.bed -b b.bed -sorted 2>&1 > /dev/null | grep "Error: Sorted input specified, " | cut -f1 -d ',' > obs
check obs exp
rm obs exp

###########################################################
#  Test that we throw an error for non-existant files
############################################################
echo -e "    intersect.new.t50...\c"
echo "Error: Unable to open file nonExistantFile.bed. Exiting." > exp
$BT intersect -a nonExistantFile.bed -b b.bed -sorted 2>&1 > /dev/null | cat - >obs
check obs exp
rm obs exp

###########################################################
#  Test that we allow empty query files
############################################################
echo -e "    intersect.new.t51...\c"
touch exp
touch dummy.txt
$BT intersect -a dummy.txt -b b.bed 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs exp dummy.txt


###########################################################
#  Test that we throw an error for unrecognized arguments
############################################################
echo -e "    intersect.new.t52...\c"
echo "***** ERROR: Unrecognized parameter: -wrongArg *****" > exp
$BT intersect -a a.bed -b b.bed -wrongArg 2>&1 > /dev/null | tail -1 > obs
check obs exp
rm obs exp


###########################################################
#  Test that we can process a Bam file with no text in 
#  it's header. 
############################################################
echo -e "    intersect.new.t53...\c"
$BT intersect -a gdc.bam -b gdc.bam -bed > obs
check obs gdc_exp
rm obs

###########################################################
#  Test that if the query is BAM, and bed output is not 
#  explicit, and they asked for an output option that is not
# valid with BAM output, an error is thrown.
############################################################
echo -e "    intersect.new.t54...\c"
echo "***** ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option. *****" > exp
$BT intersect -a a.bam -b b.bed  -wao 2>&1 > /dev/null  | tail -1 > obs
check obs exp
rm obs exp


###########################################################
# Test that gff files work correctly
############################################################
echo -e "    intersect.new.t55...\c"
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
echo -e "    intersect.new.t56...\c"
touch exp
touch dummy.txt
$BT intersect -a a.bed -b dummy.txt 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs exp dummy.txt

###########################################################
#  Test that we allow empty database files, sorted
############################################################
echo -e "    intersect.new.t57...\c"
touch exp
touch dummy.txt
$BT intersect -a a.bed -b dummy.txt -sorted 2>&1 > /dev/null | cat - > obs
check obs exp
rm obs exp dummy.txt


###########################################################
#  Test that we allow empty database files, unsorted, with
# -v (noHit) option
############################################################
echo -e "    intersect.new.t58...\c"
touch dummy.txt
$BT intersect -a a.bed -b dummy.txt -v > obs
check obs a.bed
rm obs dummy.txt


###########################################################
#  Test that an empty query with header run with -header
# option will print header
############################################################
echo -e "    intersect.new.t59...\c"
echo "#Random Header" >dummy.txt
$BT intersect -a dummy.txt -b a.bed -header > obs
check obs dummy.txt
rm obs dummy.txt


###########################################################
#  Test that an empty query with header, gzipped, that
#  runs with -header option will print header
############################################################
echo -e "    intersect.new.t60...\c"
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
echo -e "    intersect.new.t61...\c"
echo "#Random Header" | $htsutil bgzfcompress - dummy.txt.gz
echo "#Random Header" >exp
$BT intersect -a dummy.txt.gz -b a.bed -header > obs
check obs exp
rm obs dummy.txt.gz exp


###########################################################
#  Test that an empty VCF query with header that
#  runs with -header option will print header
############################################################
echo -e "    intersect.new.t62...\c"
$BT intersect -a headerOnly.vcf -b a.bed -header > obs
check obs headerOnly.vcf
rm obs


###########################################################
#  Test that files with DOS newline characters, '\r',
#  and/or extra tabs at end of line are handled
############################################################
echo -e "    intersect.new.t63...\c"
echo \
"Error: Type checker found wrong number of fields while tokenizing data line.
Perhaps you have extra TAB at the end of your line? Check with \"cat -t\""   >exp
$BT intersect -a dosLineChar_a.bed -b dosLineCharWithExtraTab_b.bed -v 2>&1 > /dev/null | cat - > obs
check exp obs
rm exp obs

###########################################################
#  Test that files with no newlines at all are handled
############################################################
echo -e "    intersect.new.t64...\c"
echo "chr17	7577068	7577157" > exp
$BT intersect -a oneRecordNoNewline.bed -b oneRecordNoNewline.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test zero length intersections in non-bam files.
############################################################
echo -e "    intersect.new.t65...\c"
echo \
"chr1	5	15	r1	chr1	9	9	m3	0
chr1	7	12	r3	chr1	9	9	m3	0" > exp
$BT intersect -a a_testZeroLen.bed -b b_testZeroLen.bed -wo > obs
check exp obs
rm exp obs

###########################################################
#  Test zero length intersections in non-bam files, -sorted
############################################################
echo -e "    intersect.new.t66...\c"
echo \
"chr1	5	15	r1	chr1	9	9	m3	0
chr1	7	12	r3	chr1	9	9	m3	0" > exp
$BT intersect -a a_testZeroLen.bed -b b_testZeroLen.bed -wo -sorted> obs
check exp obs
rm exp obs

###########################################################
#  Test vcf struct var intersection
############################################################
echo -e "    intersect.new.t67a...\c"
echo \
"19	252806	791255	G	<DEL>	70.90	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-389,-4611;END=253195;STR=+-:4;IMPRECISE;CIPOS=-2,137;CIEND=0,0;EVENT=791255;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	256900	791255	G	T	70.90	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-389,-4611;END=253195;STR=+-:4;IMPRECISE;CIPOS=-2,137;CIEND=0,0;EVENT=791255;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	260365	791256	C	<DEL>	33.71	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-680;END=261045;STR=+-:4;IMPRECISE;CIPOS=-1,257;CIEND=0,0;EVENT=791256;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=upstream_gene_variant|||ENSG00000271846|CTD-3113P16.9|ENST00000607399|||||processed_pseudogene	19	260800	791256	C	<INS>	33.71	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=680;END=261045;STR=+-:4;IMPRECISE;CIPOS=-1,257;CIEND=0,0;EVENT=791256;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=upstream_gene_variant|||ENSG00000271846|CTD-3113P16.9|ENST00000607399|||||processed_pseudogene
19	265134	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	265500	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	265986	791258	A	<DEL>	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	265500	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	265986	791258	A	<DEL>	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	266003	791258	A	C	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||" > exp
$BT intersect -a a_vcfSVtest.vcf -b b_vcfSVtest.vcf -wa -wb >obs
check exp obs
rm exp obs

echo -e "    intersect.new.t67b...\c"
echo -n "" > exp
$BT intersect -a a_vcfSVtest.vcf -b b_vcfSVtest.vcf -wa -wb -v >obs
check exp obs
rm exp obs

###########################################################
#  Test vcf struct var intersection, sorted
############################################################
echo -e "    intersect.new.t68a...\c"
echo \
"19	252806	791255	G	<DEL>	70.90	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-389,-4611;END=253195;STR=+-:4;IMPRECISE;CIPOS=-2,137;CIEND=0,0;EVENT=791255;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	256900	791255	G	T	70.90	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-389,-4611;END=253195;STR=+-:4;IMPRECISE;CIPOS=-2,137;CIEND=0,0;EVENT=791255;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	260365	791256	C	<DEL>	33.71	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-680;END=261045;STR=+-:4;IMPRECISE;CIPOS=-1,257;CIEND=0,0;EVENT=791256;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=upstream_gene_variant|||ENSG00000271846|CTD-3113P16.9|ENST00000607399|||||processed_pseudogene	19	260800	791256	C	<INS>	33.71	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=680;END=261045;STR=+-:4;IMPRECISE;CIPOS=-1,257;CIEND=0,0;EVENT=791256;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=upstream_gene_variant|||ENSG00000271846|CTD-3113P16.9|ENST00000607399|||||processed_pseudogene
19	265134	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	265500	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	265986	791258	A	<DEL>	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	265500	791257	A	<DEL>	20.25	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-558;END=265692;STR=+-:4;IMPRECISE;CIPOS=-1,196;CIEND=0,0;EVENT=791257;SUP=4;PESUP=4;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||
19	265986	791258	A	<DEL>	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||	19	266003	791258	A	C	22.15	.	TOOL=LUMPY;SVTYPE=DEL;SVLEN=-401;END=266387;STR=+-:6;IMPRECISE;CIPOS=-2,87;CIEND=0,0;EVENT=791258;SUP=6;PESUP=6;SRSUP=0;EV=PE;PRIN;CSQ=intergenic_variant||||||||||" > exp
$BT intersect -a a_vcfSVtest.vcf -b b_vcfSVtest.vcf -wa -wb -sorted >obs
check exp obs
rm exp obs

echo -e "    intersect.new.t68b...\c"
echo -n "" > exp
$BT intersect -a a_vcfSVtest.vcf -b b_vcfSVtest.vcf -wa -wb -v -sorted>obs
check exp obs
rm exp obs

###########################################################
#  Test intersect -loj with multiple databases
############################################################
echo -e "    intersect.new.t69...\c"
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
echo -e "    intersect.new.t70...\c"
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
echo -e "    intersect.new.t71...\c"
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
echo -e "    intersect.new.t72...\c"
echo \
"1	100	200	a1	ax	b	1	100	200	b1	bx	100
1	300	400	a2	ay	.	.	-1	-1	.	.	0
1	400	500	a3	az	.	.	-1	-1	.	.	0
2	500	600	a4	aq	.	.	-1	-1	.	.	0" > exp
$BT intersect -a null_a.bed -b null_b.bed null_c.bed -wao -names b c > obs
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;

###########################################################
#  Test intersect with cram without a reference provided
############################################################
echo -e "    intersect.new.t73...\c"
echo \
"FCC1MK2ACXX:2:2110:4301:28831#	99	chr1	10004	0	100M	=	10047	140	*	*	PG:Z:novoalign	AM:i:2	SM:i:2	MQ:i:0	PQ:i:219	UQ:i:0	AS:i:0	RG:Z:NCH411GBM_CD133low
FCC1MK2ACXX:2:2110:4301:28831#	147	chr1	10047	0	3S97M	=	10004	-140	*	*	PG:Z:novoalign	AM:i:2	SM:i:70	MQ:i:0	PQ:i:219	UQ:i:194	AS:i:194	RG:Z:NCH411GBM_CD133low" > exp
$BT intersect -a a.cram -b b.cram | $htsutil viewbamrecords > obs
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;


###########################################################
#  Test intersect with with pos > than 512Mb
############################################################
echo -e "    intersect.new.t74...\c"
echo \
"1	1000000004	1000000005
1	10000000004	10000000005
1	30000000000	30000000005" > exp
$BT intersect -a large_a.bed  -b large_b.bed > obs
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;

###########################################################
#  Test intersect with cram with a reference provided
############################################################
echo -e "    intersect.new.t75...\c"
echo \
"FCC1MK2ACXX:2:2110:4301:28831#	99	chr1	10004	0	100M	=	10047	140	CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT	CCCFFFFFHHHHHJJIJJJJIGIIIJIGIGJJJJIJJJJIJIIJJDIIIJJIJEHEE@GAHHGE@CEFB;A>AB=??5?B?2?<ABDAD9ABBBC9ABDA	PG:Z:novoalign	AM:i:2	SM:i:2	MQ:i:0	PQ:i:219	UQ:i:0	AS:i:0	MD:Z:100	NM:i:0	RG:Z:NCH411GBM_CD133low
FCC1MK2ACXX:2:2110:4301:28831#	147	chr1	10047	0	3S97M	=	10004	-140	AACCCTACCCCTACCCCTAACCCTACCCCTACCCCTACCCCTACCCCTACCCCTACCCCTACCCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCC	###?<55-@?250&A?882(@?795&B?8;;/?8('9'A?;8?8C;(<A@FB-/7'F?((@0D:)*9JIIHFCJJJIHFJJJIHEJJHHFFHFFFDDCCB	PG:Z:novoalign	AM:i:2	SM:i:70	MQ:i:0	PQ:i:219	UQ:i:194	AS:i:194	MD:Z:4A5A11A5A5A5A5A5A5A3A34	NM:i:10	RG:Z:NCH411GBM_CD133low">exp
CRAM_REFERENCE=test_ref.fa $BT intersect -a a.cram -b b.cram | $htsutil viewcramrecords - test_ref.fa > obs
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;

###########################################################
#  Test intersect with bam has no text header
############################################################
echo -e "    intersect.new.t76...\c"
echo \
"GA5:3:2:1710:1301#0	0	dummy_chr	1279	18	76M	*	0	0	GAACTCCTGACCTCAGGTGATCTGCCCGCCTTGGCCTCCCAAAGTGCTGGAATTACAGGCATGAGCCACCGTGCCC	6->;B==?B?AAA??9B<AA?@A==AA9A?<A?<&:A?=<9<?>19;?A=9:?999;=4=6A9/8;4==;1';;3=	XT:A:U	NM:i:0	X0:i:1	X1:i:3	XM:i:0	XO:i:0	XG:i:0	MD:Z:76	XA:Z:chr14,+56265423,76M,1;chr14,+91530561,76M,1;chr2,+230545178,76M,1;
GA5:3:24:462:583#0	0	dummy_chr	128882	37	76M	*	0	0	TAAAAAAAGGACAGTGACGCACCTTGTATAGCGATGTGTCATCTAAAACATCTATTCAAAGAACAGAAGACTCACC	BBABBBB@BBBBBA?BBBBBBAA@AB?BBBA@B=>@?@;@?@?A>=1=?>??><???;;>?;;;?>8999;;9;;<	XT:A:U	NM:i:3	X0:i:1	X1:i:0	XM:i:3	XO:i:0	XG:i:0	MD:Z:6C21G30G16
GA5:3:29:1241:1653#0	0	dummy_chr	5591013	37	76M	*	0	0	TCATGCACACACAGACAGCTGTCGGGGGATGCATGCCAACCAGAGGGGCCACACATATACCGTGTTGATGGGACAG	BBBBBBACBCBBBBBBBBBBBABBABBB@AA>@?@?@>=????<???:?915399<=5=<==5=559545432353	XT:A:U	NM:i:1	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:58T17
GA5:3:33:1591:303#0	0	dummy_chr	11880048	37	76M	*	0	0	CTCGCCTGGGCCCGGTAAAGCCCCCACGTAGCCCCAGCCAGCCTGGAACATGCTTCCTGAGCTCCCAGCTCTTGGT	BBBBBBBBBBB@BBB?BBBA@A>B@B?AABB?@>?BA?>AB;??A?6;=AAAAA=>3=9;@;===6,=?;;5==;6	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:76
GA5:3:31:677:1537#0	16	dummy_chr	11880931	37	76M	*	0	0	GAGGGTTTGAGAGAGCAGCCAGGAGAGCTTAGGGTCTCAGGGTGTCCCAGACCCCGACACCGGCCAGTGGCGGAAG	###################>9==9=9?>????>>?>8>>??AA?>A??>A4AA?AA?B=ABBBBBBBBBBBBBBBB	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:76
GA5:3:49:1480:1116#0	16	dummy_chr	11913868	37	76M	*	0	0	GGGAGGAGGCCAGGACTTCAGGGACCCACAGCCATCACCTCCCTCCCCTGCCCCCTACACACCAACTCTCTGGAAA	#################################44:4=944==;=???>=?>==??=A=A;ABA?A?AABAAABBB	XT:A:U	NM:i:1	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0T75
GA5:3:61:213:1812#0	16	dummy_chr	13030396	37	76M	*	0	0	GGTCCGGCGGGGTCGGACTGGACCAGCTGTTGGGCTTTGTTTGCTCTTTTTACGAATTGAAAAACTGAAGCCAGGA	/=81,5948=485=4,),1;;7:87:6=;;@@AB=C8A@@BAB=>5>BBBB>BBAAA9ABA@B4B;BBBBBBBBCB	XT:A:U	NM:i:1	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0T75
GA5:3:116:1581:552#0	16	dummy_chr	15055984	37	76M	*	0	0	AGAAAGCCTAAGGTCAGGGTGCCAGCAGGTTTGGTGTCTGGTGAGGTACCCATCTCTGCTTCTAAGGCAGAGCCTT	48887429,3=;98<8<8@;<=?8@??98@@@<=AA>@@?B?A@@BA6BA@=@BABBB???B@BBBBBABBCBB?B	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:76" > exp
$BT intersect -a notexthdr.bam -b notexthdr.bam | $htsutil viewbamrecords > obs
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;

###########################################################
#  Test intersect with 2G io buffer, shouldn't overflow
############################################################
echo -e "    intersect.new.t77...\c"
echo -n "" > exp
$BT intersect -iobuf 2G -ubam -S -u -sorted -b a.bam -a a.bed >obs
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;



###########################################################
#  Test intersect preserve the text header in bam
############################################################
echo -e "    intersect.new.t78...\c"
echo -e "@HD	VN:1.5	SO:coordinate" > exp
echo "@HD	VN:1.5	SO:coordinate" | $htsutil samtobam - - | $BT intersect -a /dev/stdin -b b.bed | $htsutil viewbam >obs
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;


###########################################################
#  Test intersect recognize "*" as a valid strand char
############################################################
echo -e "    intersect.new.t79...\c"
echo	-e	"GL000008.2	118286	119434	Pituitary_00000002	999	*
GL000008.2	158006	158759	Pituitary_00000012	999	+
GL000009.2	78532	79175	Pituitary_00000033	999	*"	>	exp
$BT  intersect -s -v -a strand_with_star.a.bed -b strand_with_star.b.bed > obs 
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;


###########################################################
#  Test utral bed4 intersection
############################################################
echo -e "    intersect.new.t80...\c"
cp ultra-long-bed4.bed exp
$BT intersect -a ultra-long-bed4.bed -b ultra-long-bed4.bed > obs
check exp obs
rm exp obs
[[ $FAILURES -eq 0 ]] || exit 1;

###########################################################
#  Test -C that got -b parameter first
############################################################
echo -e "    intersect.new.t81...\c"
echo	-e	"chr1	10	20	a1	1	+	0
chr1	100	200	a2	2	-	2"	>	exp
$BT intersect -b b.bed -a a.bed -C > obs
check exp obs
rm exp obs
