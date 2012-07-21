###############
5.12 maskFastaFromBed
###############
**maskFastaFromBed** masks sequences in a FASTA file based on intervals defined in a feature file. The
headers in the input FASTA file must exactly match the chromosome column in the feature file. This
may be useful fro creating your own masked genome file based on custom annotations or for masking all
but your target regions when aligning sequence data from a targeted capture experiment.


==========================================================================
5.12.1 Usage and option summary
==========================================================================
Usage:
::
  maskFastaFromBed [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>
  
NOTE: The input and output FASTA files must be different.

===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-soft**				         Soft-mask (that is, convert to lower-case bases) the FASTA sequence. *By default, hard-masking (that is, conversion to Ns) is performed*.							                   
===========================      ===============================================================================================================================================================================================================






==========================================================================
5.12.2 Default behavior
==========================================================================
**maskFastaFromBed** will mask a FASTA file based on the intervals in a BED file. The newly masked
FASTA file is written to the output FASTA file.

For example:
::
  cat test.fa
  >chr1
  AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

  cat test.bed
  chr1 5 10

  maskFastaFromBed -fi test.fa -bed test.bed -fo test.fa.out
  
  cat test.fa.out
  >chr1
  AAAAANNNNNCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG


==========================================================================
5.12.3 Soft-masking the FASTA file.
==========================================================================
Using the **-soft** option, one can optionally "soft-mask" the FASTA file.

For example:
::
  cat test.fa
  >chr1
  AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

  cat test.bed
  chr1 5 10

  maskFastaFromBed -fi test.fa -bed test.bed -fo test.fa.out -soft

  cat test.fa.out
  >chr1
  AAAAAaaaccCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG
