###############
5.4 bamToBed
###############

**bamToBed** is a general purpose tool that will convert sequence alignments in BAM format to either
BED6, BED12 or BEDPE format. This enables one to convert BAM files for use with all of the other
BEDTools. The CIGAR string is used to compute the alignment end coordinate in an "ungapped"
fashion. That is, match ("M"), deletion ("D"), and splice ("N") operations are observed when computing
alignment ends.

============================================
5.4.1 Usage and option summary
============================================
**Usage:**
::
  bamToBed [OPTIONS] -i <BAM>
  

======================           =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
Option                              Description
======================           =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
**-bedpe**				         Write BAM alignments in BEDPE format. Only one alignment from paired-end reads will be reported. Specifically, it each mate is aligned to the same chromosome, the BAM alignment reported will be the one where the BAM insert size is greater than zero. When the mate alignments are interchromosomal, the lexicographically lower chromosome will be reported first. Lastly, when an end is unmapped, the chromosome and strand will be set to "." and the start and end coordinates will be set to -1. *By default, this is disabled and the output will be reported in BED format*.								 								 
								 **NOTE: When using this option, it is required that the BAM file is sorted/grouped by the read name. This allows bamToBed to extract correct alignment coordinates for each end based on their respective CIGAR strings. It also assumes that the alignments for a given pair come in groups of twos. There is not yet a standard method for reporting multiple alignments using BAM. bamToBed will fail if an aligner does not report alignments in pairs**.							 
                                 BAM files may be piped to bamToBed by specifying "-i stdin". See example below.
**-bed12**					     Write "blocked" BED (a.k.a. BED12) format. This will convert "spliced" BAM alignments (denoted by the "N" CIGAR operation) to BED12.
**-ed**					         Use the "edit distance" tag (NM) for the BED score field. Default for BED is to use mapping quality. Default for BEDPE is to use the *minimum* of the two mapping qualities for the pair. When -ed is used with -bedpe, the total edit distance from the two mates is reported.                                            
**-tag**					     Use other *numeric* BAM alignment tag for BED score. Default for BED is to use mapping quality. Disallowed with BEDPE output.
**-color**					     An R,G,B string for the color used with BED12 format. Default is (255,0,0).                              
**-split**					     Report each portion of a "split" BAM (i.e., having an "N" CIGAR operation) alignment as a distinct BED intervals.			            
======================           =========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

By default, each alignment in the BAM file is converted to a 6 column BED. The BED "name" field is
comprised of the RNAME field in the BAM alignment. If mate information is available, the mate (e.g.,
"/1" or "/2") field will be appended to the name. The "score" field is the mapping quality score from the
BAM alignment, unless the **-ed** option is used.

Examples:
::
  bamToBed -i reads.bam | head -5
  chr7   118970079   118970129   TUPAC_0001:3:1:0:1452#0/1   37   -
  chr7   118965072   118965122   TUPAC_0001:3:1:0:1452#0/2   37   +
  chr11  46769934    46769984    TUPAC_0001:3:1:0:1472#0/1   37   -
  
  bamToBed -i reads.bam -tag NM | head -5
  chr7   118970079   118970129   TUPAC_0001:3:1:0:1452#0/1   1    -
  chr7   118965072   118965122   TUPAC_0001:3:1:0:1452#0/2   3    +
  chr11  46769934    46769984    TUPAC_0001:3:1:0:1472#0/1   1    -
  
  bamToBed -i reads.bam -bedpe | head -3
  chr7   118965072   118965122   chr7   118970079   118970129
         TUPAC_0001:3:1:0:1452#0 37     +     -
  chr11  46765606    46765656    chr11  46769934    46769984
         TUPAC_0001:3:1:0:1472#0 37     +     -
  chr20  54704674    54704724    chr20  54708987    54709037
         TUPAC_0001:3:1:1:1833#0 37     +    

		 
One can easily use samtools and bamToBed together as part of a UNIX pipe. In this example, we will
only convert properly-paired (BAM flag == 0x2) reads to BED format.
::
  samtools view -bf 0x2 reads.bam | bamToBed -i stdin | head
  chr7   118970079   118970129   TUPAC_0001:3:1:0:1452#0/1   37   -
  chr7   118965072   118965122   TUPAC_0001:3:1:0:1452#0/2   37   +
  chr11  46769934    46769984    TUPAC_0001:3:1:0:1472#0/1   37   -
  chr11  46765606    46765656    TUPAC_0001:3:1:0:1472#0/2   37   +
  chr20  54704674    54704724    TUPAC_0001:3:1:1:1833#0/1   37   +
  chr20  54708987    54709037    TUPAC_0001:3:1:1:1833#0/2   37   -
  chrX   9380413     9380463     TUPAC_0001:3:1:1:285#0/1    0    -
  chrX   9375861     9375911     TUPAC_0001:3:1:1:285#0/2    0    +
  chrX   131756978   131757028   TUPAC_0001:3:1:2:523#0/1    37   +
  chrX   131761790   131761840   TUPAC_0001:3:1:2:523#0/2    37   -

  
==================================================================
5.4.2 (-split)Creating BED12 features from "spliced" BAM entries. 
==================================================================
bamToBed will, by default, create a BED6 feature that represents the entire span of a spliced/split
BAM alignment. However, when using the **-split** command, a BED12 feature is reported where BED
blocks will be created for each aligned portion of the sequencing read.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             
  Exons       ***************                                    **********
  
  BED/BAM A      ^^^^^^^^^^^^....................................^^^^
  
  Result      ===============                                    ====
  
