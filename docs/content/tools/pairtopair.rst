.. _pairtopair:

###############
*pairtopair*
###############
**pairToPair** compares two BEDPE files in search of overlaps where each end of a BEDPE feature in A
overlaps with the ends of a feature in B. For example, using pairToPair, one could screen for the exact
same discordant paired-end alignment in two files. This could suggest (among other things) that the
discordant pair suggests the same structural variation in each file/sample.


================================
Usage and option summary
================================
**Usage:**

::

  pairToPair [OPTIONS] -a <BEDPE> -b <BEDPE>
  
  
===========================      =========================================================================================================================================================
Option                           Description
===========================      =========================================================================================================================================================
**-a**				             BEDPE file A. Each feature in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe.
**-b**					         BEDPE file B. Use "stdin" if passing B with a UNIX pipe.
**-f** 				             Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
**-is** 				         Force "strandedness". That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
**-type**					     Approach to reporting overlaps between BEDPE and BED.
                                 | **either** Report overlaps if either ends of A overlap B.	
								     
								 
								 | **neither** Report A if neither end of A overlaps B.
							
								 
								 | **both** Report overlaps if both ends of A overlap B.   -*Default behavior.*
===========================      =========================================================================================================================================================





================================
Default behavior
================================
By default, a BEDPE feature from A will be reported if *both* ends overlap a feature in the BEDPE B
file. If strand information is present for the two BEDPE files, it will be further required that the
overlaps on each end be on the same strand. This way, an otherwise overlapping (in terms of genomic
locations) F/R alignment will not be matched with a R/R alignment.

Default: Report A if *both* ends overlaps B.

::

  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^^^^^                                          ^^^^^^
  
  Result              =====.................................=====


Default when strand information is present in both BEDPE files: Report A if *both* ends overlaps B *on
the same strands*.

::

  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE A         >>>>>.................................>>>>>
  
  BEDPE B            <<<<<.............................>>>>>
  
  Result
  
  
  
  BEDPE A         >>>>>.................................>>>>>
  
  BEDPE B            >>>>>.............................>>>>>
  
  Result          >>>>>.................................>>>>> 


  
==================================================
``-type neither`` Optional overlap requirements 
==================================================
Using then **-type neither, pairToPair** will only report A if *neither* end overlaps with a BEDPE
feature in B.

**-type neither**: Report A only if *neither* end overlaps B.

::

  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B         ^^^^^^^^......................................^^^^^^
  
  Result             
  
  
  
  BEDPE/BAM A         *****.................................*****
  
  BED File B    ^^^^................................................^^^^^^
  
  Result              =====.................................=====
  
  
