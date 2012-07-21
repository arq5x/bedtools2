###############
5.2 pairToBed
###############
**pairToBed** compares each end of a BEDPE feature or a paired-end BAM alignment to a feature file in
search of overlaps.

**NOTE: pairToBed requires that the BAM file is sorted/grouped by the read name. This
allows pairToBed to extract correct alignment coordinates for each end based on their
respective CIGAR strings. It also assumes that the alignments for a given pair come in
groups of twos. There is not yet a standard method for reporting multiple alignments
using BAM. pairToBed will fail if an aligner does not report alignments in pairs.**

==========================================================================
5.2.1 Usage and option summary
==========================================================================
**Usage:**
::
  pairToBed [OPTIONS] [-a <BEDPE> || -abam <BAM>] -b <BED/GFF/VCF>
  
  
===========================      =========================================================================================================================================================
Option                           Description
===========================      =========================================================================================================================================================
**-a**				             BEDPE file A. Each feature in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe. Output will be in BEDPE format.
**-b**					         BED file B. Use "stdin" if passing B with a UNIX pipe.
**-abam**					     BAM file A. Each end of each BAM alignment in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe: For example: samtools view 每b <BAM> | pairToBed 每abam stdin 每b genes.bed | samtools view -                                                  
**-ubam**					     Write uncompressed BAM output. The default is write compressed BAM output.
**-bedpe**					     When using BAM input (-abam), write output as BEDPE. The default is to write output in BAM when using -abam. For example: pairToBed 每abam reads.bam 每b genes.bed 每bedpe                              
**-ed**					         Use BAM total edit distance (NM tag) for BEDPE score. Default for BEDPE is to use the *minimum* of the two mapping qualities for the pair. When -ed is used the *total* edit distance from the two mates is reported as the score.
**-f** 				             Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
**-s** 				             Force "strandedness". That is, only report hits in B that overlap A on the **same** strand. By default, overlaps are reported without respect to strand.
**-type**					     
                                 Approach to reporting overlaps between BEDPE and BED.

                                 
								 **either-** Report overlaps if either end of A overlaps B.	
								 
								 - *Default*
								 
								 **neither-** Report A if neither end of A overlaps B.
								 
								 **xor-** Report overlaps if one and only one end of A overlaps B.
								 
								 **both-** Report overlaps if both ends of A overlap B.
								 
								 **notboth-** Report overlaps if neither end or one and only one end of A overlap B.
								 
								 **ispan-** Report overlaps between [end1, start2] of A and B.	
								 
								 - Note: If chrom1 <> chrom2, entry is ignored.
								  
							     **ospan-** Report overlaps between [start1, end2] of A and B.
								 
								 - Note: If chrom1 <> chrom2, entry is ignored.
									   
								 **notispan-**  Report A if ispan of A doesn't overlap B.
								 - Note: If chrom1 <> chrom2, entry is ignored.
												
								 **notospan-**  Report A if ospan of A doesn't overlap B.
								 - Note: If chrom1 <> chrom2, entry is ignored.
===========================      =========================================================================================================================================================



==========================================================================
5.2.2 Default behavior
==========================================================================
By default, a BEDPE / BAM feature will be reported if *either* end overlaps a feature in the BED file.
In the example below, the left end of the pair overlaps B yet the right end does not. Thus, BEDPE/
BAM A is reported since the default is to report A if either end overlaps B.

Default: Report A if *either* end overlaps B.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^^^^^                                          ^^^^^^
  
  Result              =====.................................=====

  
==========================================================================
5.2.3 (-type)Optional overlap requirements 
==========================================================================
Using then **-type** option, **pairToBed** provides several other overlap requirements for controlling how
overlaps between BEDPE/BAM A and BED B are reported. The examples below illustrate how each
option behaves.

**-type both**: Report A only if *both* ends overlap B.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^^^^^                                          ^^^^^^
  
  Result
  
  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^^^^^                                   ^^^^^^
  
  Result              =====.................................=====
  
  
**-type neither**: Report A only if *neither* end overlaps B.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^^^^^                                          ^^^^^^
  
  Result
  
  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B   ^^^^                                                  ^^^^^^
  
  Result              =====.................................=====
  
  
**-type xor**: Report A only if *one and only one* end overlaps B.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^^^^^                                          ^^^^^^
  
  Result              =====.................................=====      

  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^                                   ^^^^^^
  
  Result  
  
  
**-type notboth**: Report A only if *neither end* **or** *one and only one* end overlaps B. Thus "notboth"
includes what would be reported by "neither" and by "xor".
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^^^^^                                          ^^^^^^
  
  Result              =====.................................=====     
  
  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B     ^^^                                               ^^^^^^
  
  Result              =====.................................=====   
  
  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^                                   ^^^^^^
  
  Result  
  
  
**-type ispan**: Report A if it's "*inner span*" overlaps B. Applicable only to intra-chromosomal features.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
                Inner span |-------------------------------|
				
  BEDPE/BAM A         *****.................................*****
  
  BED File B                         ^^^^^^^^                       
  
  Result              =====.................................=====    
  
  
  
  BEDPE/BAM A         =====.................................=====
  
  BED File B         ====
  
  Result
  

**-type ospan**: Report A if it's "*outer span*" overlaps B. Applicable only to intra-chromosomal features.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
          Outer span  |-----------------------------------------|
		  
  BEDPE/BAM A         *****.................................*****
  
  BED File B             ^^^^^^^^^^^^
  
  Result              =====.................................=====    
  
  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B     ^^^^
  
  Result
  
  
**-type notispan**: Report A only if it's "*inner span*" does not overlap B. Applicable only to intrachromosomal
features.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
                Inner span |-------------------------------|
				
  BEDPE/BAM A         *****.................................*****
  
  BED File B                         ^^^^^^^^
  
  Result             
  
  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^
  
  Result              =====.................................=====
  
  
**-type notospan**: Report A if it's "*outer span*" overlaps B. Applicable only to intra-chromosomal
features.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
          Outer span  |-----------------------------------------|
		  
  BEDPE/BAM A         *****.................................*****
                         
  BED File B             ^^^^^^^^^^^^
  
  Result                 


  
  BEDPE/BAM A         *****.................................*****
  
  BED File B     ^^^^
  
  Result              =====.................................===== 
  
  

==========================================================================
5.2.4 (-f)Requiring a minimum overlap fraction 
==========================================================================
By default, **pairToBed** will report an overlap between A and B so long as there is at least one base
pair is overlapping on either end. Yet sometimes you may want to restrict reported overlaps between A
and B to cases where the feature in B overlaps at least X% (e.g. 50%) of A. The **每f** option does exactly
this. The **-f** option may also be combined with the -type option for additional control. For example,
combining **-f 0.50** with **-type both** requires that both ends of A have at least 50% overlap with a
feature in B.

For example, report A only at least 50% of one of the two ends is overlapped by B.
:: 
  pairToBed -a A.bedpe -b B.bed -f 0.5


  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^                                           ^^^^^^
  
  Result                  

  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^                                         ^^^^^^
  
  Result              =====.................................=====
  

  
==========================================================================
5.2.5 (-s)Enforcing "strandedness" 
==========================================================================
By default, **pairToBed** will report overlaps between features even if the features are on opposing
strands. However, if strand information is present in both files and the **"-s"** option is used, overlaps will
only be reported when features are on the same strand.

For example, report A only at least 50% of one of the two ends is overlapped by B.
::
  pairToBed -a A.bedpe -b B.bed -s
  

  
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         >>>>>.................................<<<<<
  
  BED File B         <<                                           >>>>>
  
  Result
  
  
  
  BEDPE/BAM A         >>>>>.................................<<<<<
  
  BED File B         >>                                          >>>>>
  
  Result              >>>>>.................................<<<<<
  
  
==========================================================================
5.2.6 (-abam)Default is to write BAM output when using BAM input 
==========================================================================
When comparing *paired* alignments in BAM format (**-abam**) to features in BED format (**-b**),
**pairToBed** will , by default, write the output in BAM format. That is, each alignment in the BAM
file that meets the user's criteria will be written (to standard output) in BAM format. This serves as a
mechanism to create subsets of BAM alignments are of biological interest, etc. Note that both
alignments for each aligned pair will be written to the BAM output.

For example:
::
  pairToBed 每abam pairedReads.bam 每b simreps.bed | samtools view - | head -4

  JOBU_0001:3:1:4:1060#0 99 chr10 42387928 29 50M = 42393091 5 2 1 3
  AA A A A C G G A A T T A T C G A A T G G A A T C G A A G A G A A T C T T C G A A C G G A C C C G A
  dcgggggfbgfgdgggggggfdfgggcggggfcggcggggggagfgbggc XT:A:R NM:i:5 SM:i:0 AM:i:0 X0:i:3 X 1 : i :
  3 XM:i:5 XO:i:0 XG:i:0 MD:Z:0T0C33A5T4T3
  JOBU_0001:3:1:4:1060#0 147 chr10 42393091 0 50M = 42387928 - 5 2 1 3
  AAATGGAATCGAATGGAATCAACATCAAATGGAATCAAATGGAATCATTG K g d c g g d e c d g
  \d`ggfcgcggffcgggc^cgfgccgggfc^gcdgg\bg XT:A:R NM:i:2 SM:i:0 AM:i:0 X0:i:3 X1:i:13 XM:i:2 X O : i :
  0 XG:i:0 MD:Z:21T14G13
  JOBU_0001:3:1:8:446#0 99 chr10 42388091 9 50M = 42392738 4 6 9 7
  GAATCGACTGGAATCATCATCGGATGGAAATGAATGGAATAATCATCGAA f _ O f f ` ] I e Y f f ` f f e d d c f e f c P ` c _ W \ \ R _ ]
  _BBBBBBBBBBBBBBBB XT:A:U NM:i:4 SM:i:0 AM:i:0 X0:i:1 X1:i:3 XM:i:4 XO:i:0 XG:i:0 M D : Z :
  7A22C9C2T6
  JOBU_0001:3:1:8:446#0 147 chr10 42392738 9 50M = 42388091 - 4 6 9 7
  TTATCGAATGCAATCGAATGGAATTATCGAATGCAATCGAATAGAATCAT df^ffec_JW[`MWceRec``fee`dcecfeeZae`c]
  f^cNeecfccf^ XT:A:R NM:i:1 SM:i:0 AM:i:0 X0:i:2 X1:i:2 XM:i:1 XO:i:0 XG:i:0 MD:Z:38A11
  
  
  
==========================================================================
5.2.7 (-bedpe)Output BEDPE format when using BAM input 
==========================================================================
When comparing *paired* alignments in BAM format (**-abam**) to features in BED format (**-b**),
**pairToBed** will optionally write the output in BEDPE format. That is, each alignment in the BAM
file is converted to a 10 column BEDPE feature and if overlaps are found (or not) based on the user's
criteria, the BAM alignment will be reported in BEDPE format. The BEDPE "name" field is comprised
of the RNAME field in the BAM alignment. The "score" field is the mapping quality score from the
BAM alignment.

For example:
::
  pairToBed 每abam pairedReads.bam 每b simreps.bed -bedpe | head -5
  chr10 42387927     42387977    chr10   42393090   42393140
        JOBU_0001:3:1:4:1060#0   29      +     -
  chr10 42388090 42388140        chr10   42392737   42392787
        JOBU_0001:3:1:8:446#0    9       +     -
  chr10 42390552 42390602        chr10   42396045   42396095
        JOBU_0001:3:1:10:1865#0  9       +     -
  chrX  139153741 139153791      chrX    139159018  139159068
        JOBU_0001:3:1:14:225#0   37      +     -
  chr4  9236903 9236953          chr4    9242032    9242082
        JOBU_0001:3:1:15:1362#0  0       +     -
