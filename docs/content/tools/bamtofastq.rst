.. _bamtofastq:

###############
*bamtofastq*
###############
``bedtools bamtofastq`` is a conversion utility for extracting FASTQ records
from sequence alignments in BAM format. 


.. note::

    If you are using CRAM as input, you will need to specify
    the *full path* describing the location of the relevant reference genome in FASTA format via the CRAM_REFERENCE environment variable. For example:

        `export CRAM_REFERENCE=/path/to/ref/g1k_v37_decoy.fa`


==========================================================================
Usage and option summary
==========================================================================
**Usage**:
::

    bedtools bamtofastq [OPTIONS] -i <BAM> -fq <FASTQ>

**(or)**:
::

    bamToFastq [OPTIONS] -i <BAM> -fq <FASTQ>



.. tabularcolumns:: |p{4.5cm}|p{8.5cm}|

=============   ================================================================
Option          Description
=============   ================================================================
**-fq2**        FASTQ for second end.  Used if BAM contains paired-end data.
                BAM should be sorted by query name 
                (``samtools sort -n -o aln.qsort.bam aln.bam``) if creating 
                paired FASTQ with this option.
**-tags**       Create FASTQ based on the mate info in the BAM R2 and Q2 tags.
=============   ================================================================


==========================================================================
Default behavior
==========================================================================
By default, each alignment in the BAM file is converted to a FASTQ record
in the ``-fq`` file. The order of the records in the resulting FASTQ exactly
follows the order of the records in the BAM input file.

.. code-block:: bash

  $ bedtools bamtofastq -i NA18152.bam -fq NA18152.fq
  
  $ head -8 NA18152.fq
  @NA18152-SRR007381.35051
  GGAGACATATCATATAAGTAATGCTAGGGTGAGTGGTAGGAAGTTTTTTCATAGGAGGTGTATGAGTTGGTCGTAGCGGAATCGGGGGTATGCTGTTCGAATTCATAAGAACAGGGAGGTTAGAAGTAGGGTCTTGGTGACAAAATATGTTGTATAGAGTTCAGGGGAGAGTGCGTCATATGTTGTTCCTAGGAAGATTGTAGTGGTGAGGGTGTTTATTATAATAATGTTTGTGTATTCGGCTATGAAGAATAGGGCGAAGGGGCCTGCGGCGTATTCGATGTTGAAGCCTGAGACTAGTTCGGACTCCCCTTCGGCAAGGTCGAA
  +
  <<<;;<;<;;<;;;;;;;;;;;;<<<:;;;;;;;;;;;;;;;;::::::;;;;<<;;;;;;;;;;;;;;;;;;;;;;;;;;;;<<<<<;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;<<;;;;;:;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;<<<;;;;;;;;;;<<<<<<<<;;;;;;;;;:;;;;;;;;;;;;;;;;;;;:;;;;8;;8888;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;8966689666666299866669:899
  @NA18152-SRR007381.637219
  AATGCTAGGGTGAGTGGTAGGAAGTTTTTTCATAGGAGGTGTATGAGTTGGTCGTAGCGGAATCGGGGGTATGCTGTTCGAATTCATAAGAACAGGGAGGTTAGAAGTAGGGTCTTGGTGACAAAATATGTTGTATAGAGTTCAGGGGAGAGTGCGTCATATGTTGTTCCTAGGAAGATTGTAGTGGTGAGGGTGTTTATTATAATAATGTTTGTGTATTCGGCTATGAAGAATAGGGCGAAGGGGCCTGCGGCGTATTCGATGTTGAAGCCTGAGACTAGTTCGGACTCCCCTTCCGGCAAGGTCGAA
  +
  <<<<<<<<<<;;<;<;;;;<<;<888888899<;;;;;;<;;;;;;;;;;;;;;;;;;;;;;;;<<<<<;;;;;;;;;<;<<<<<;;;;;;;;;;;;;<<<<;;;;;;;:::;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;<<<<;;;;;;;;;;;;;;;;;;;;;;;<;;;;;;;;;;;;;;;;;;;;;;<888<;<<;;;;<<<<<<;;;;;<<<<<<<<;;;;;;;;;:;;;;888888899:::;;8;;;;;;;;;;;;;;;;;;;99;;99666896666966666600;96666669966



==========================================================================
``-fq2`` Creating two FASTQ files for paired-end sequences.
==========================================================================
If your BAM alignments are from paired-end sequence data, one can use the
``-fq2`` option to create two distinct FASTQ output files --- one for 
end 1 and one for end 2.

.. note::

    When using this option, it is required that the BAM 
    file is sorted/grouped by the read name. This keeps the resulting records
    in the two output FASTQ files in the same order. One can sort the BAM
    file by query name with ``samtools sort -n -o aln.qsort.bam aln.bam``.


.. code-block:: bash

  $ samtools sort -n -o aln.qsort.bam aln.bam
  
  $ bedtools bamtofastq -i aln.qsort.bam \
                        -fq aln.end1.fq \
                        -fq2 aln.end2.fq
                        
  $ head -8 aln.end1.fq
  @SRR069529.2276/1
  CAGGGAGAAGGAGGTAGGAAAGAGAAAGGACCAGGGAGGGGCGCATACACAGGACGCTCCGTGCGGTGATAGCAGCACCACACTGTGTTCAGTCGTCTGGC
  +
  =;@>==###############################################################################################
  @SRR069529.2406/1
  GCTGGGAAAAGGATTCAGGATGTTGGTTTCTATCTTTGAGTTGCTGCTGTGCGGCTGTCCCTACACTCGCAGTACCCCTCGGACACCGTCTACTGTGGAGG
  +
  =5@><<:?<?
  
  $ head -8 aln.end2.fq
  @SRR069529.2276/2
  AGACCCAGAGAGGGACAGGATCTGTCCCAGATCATAAAATAGGGGGAGTGCTCCGTAGAGGCGTGCGCGGTGGCACCGTGCAGTAGTACGGGTGAGCGGGG
  +
  #####################################################################################################
  @SRR069529.2406/2
  TTCCCTACCCCTGGGGTCAGGGACTACAGCCAAGGGGAGAACTTTAGCAAGTAGACGTTAGTTATTTTGATTCCAGTGGGGACGCGCGTGTAGCGAGTTGT
  +
  @>=AABB?AAACABBA>@?AAAA>B@@AB@AA:B@AA@??#############################################################


