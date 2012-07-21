###############
5.11 fastaFromBed
###############
**fastaFromBed** extracts sequences from a FASTA file for each of the intervals defined in a BED file.
The headers in the input FASTA file must exactly match the chromosome column in the BED file.

==========================================================================
5.11.1 Usage and option summary 
==========================================================================
Usage:
::
  fastaFromBed [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>

===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-name**				         Use the "name" column in the BED file for the FASTA headers in the output FASTA file.								 
**-tab**					     Report extract sequences in a tab-delimited format instead of in FASTA format.
**-s**                           Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. *Default: strand information is ignored*.
===========================      ===============================================================================================================================================================================================================







==========================================================================
5.11.2 Default behavior
==========================================================================
**fastaFromBed** will extract the sequence defined by the coordinates in a BED interval and create a
new FASTA entry in the output file for each extracted sequence. By default, the FASTA header for each
extracted sequence will be formatted as follows: "<chrom>:<start>-<end>".

For example:
::
  $ cat test.fa
  >chr1
  AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

  cat test.bed
  chr1 5 10

  fastaFromBed -fi test.fa -bed test.bed -fo test.fa.out

  cat test.fa.out
  >chr1:5-10
  AAACC



  
==========================================================================
5.11.3 Using the BED "name" column as a FASTA header.
==========================================================================
Using the **-name** option, one can set the FASTA header for each extracted sequence to be the "name"
columns from the BED feature.

For example:
::
  cat test.fa
  >chr1
  AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

  cat test.bed
  chr1 5 10 myseq

  fastaFromBed -fi test.fa -bed test.bed -fo test.fa.out -name

  cat test.fa.out
  >myseq
  AAACC










==========================================================================
5.11.4 Creating a tab-delimited output file in lieu of FASTA output.
==========================================================================
Using the **-tab** option, the **-fo** output file will be tab-delimited instead of in FASTA format.

For example:
::
  cat test.fa
  >chr1
  AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

  cat test.bed
  chr1 5 10 myseq

  fastaFromBed -fi test.fa -bed test.bed -fo test.fa.out.tab -name -tab

  cat test.fa.out
  myseq AAACC
  
  
  
==========================================================================
5.11.5 (-s)Forcing the extracted sequence to reflect the requested strand 
==========================================================================
**fastaFromBed** will extract the sequence in the orientation defined in the strand column when the "-s"
option is used.

For example:
::
  cat test.fa
  >chr1
  AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

  cat test.bed
  chr1 20 25 forward 1 +
  chr1 20 25 reverse 1 -

  fastaFromBed -fi test.fa -bed test.bed -s -name -fo test.fa.out

  cat test.fa.out
  >forward
  CGCTA
  >reverse
  TAGCG
