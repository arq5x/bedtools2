.. _annotate:

###############
*annotate*
###############
``bedtools annotate``, well, annotates one BED/VCF/GFF file with the coverage 
and number of overlaps observed from multiple other BED/VCF/GFF files. 
In this way, it allows one to ask to what degree one feature coincides with 
multiple other feature types with a single command.

==========================================================================
Usage and option summary
==========================================================================
**Usage**:
::

  bedtools annotate [OPTIONS] -i <BED/GFF/VCF> -files FILE1 FILE2 FILE3 ... FILEn

**(or)**:
::

  annotateBed [OPTIONS] -i <BED/GFF/VCF> -files FILE1 FILE2 FILE3 ... FILEn
  
  
===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-names**				         A list of names (one per file) to describe each file in -i. These names will be printed as a header line. 
**-counts**					     Report the count of features in each file that overlap -i. Default behavior is to report the fraction of -i covered by each file.
**-both**                        Report the count of features followed by the % coverage for each annotation file. Default is to report solely the fraction of -i covered by each file.
**-s**                           Force strandedness. That is, only include hits in A that overlap B on the same strand. By default, hits are included without respect to strand.
**-S**	                         Require different strandedness.  That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.
===========================      ===============================================================================================================================================================================================================

==========================================================================
Default behavior - annotate one file with coverage from others.
==========================================================================
By default, the fraction of each feature covered by each annotation file is 
reported after the complete feature in the file to be annotated.

.. code-block:: bash

  $ cat variants.bed
  chr1 100  200   nasty 1  -
  chr2 500  1000  ugly  2  +
  chr3 1000 5000  big   3  -

  $ cat genes.bed
  chr1 150  200   geneA 1  +
  chr1 175  250   geneB 2  +
  chr3 0    10000 geneC 3  -

  $ cat conserve.bed
  chr1 0    10000 cons1 1  +
  chr2 700  10000 cons2 2  -
  chr3 4000 10000 cons3 3  +

  $ cat known_var.bed
  chr1 0    120   known1   -
  chr1 150  160   known2   -
  chr2 0    10000 known3   +

  $ bedtools annotate -i variants.bed -files genes.bed conserve.bed known_var.bed
  chr1	100	200	nasty	1	-	0.500000	1.000000	0.300000	
  chr2	500	1000	ugly	2	+	0.000000	0.600000	1.000000	
  chr3	1000	5000	big	3	-	1.000000	0.250000	0.000000


==========================================================================
``-count`` Report the count of hits from the annotation files
==========================================================================

.. code-block:: bash

  $ bedtools annotate -counts -i variants.bed -files genes.bed conserve.bed known_var.bed
  chr1	100	200	nasty	1	-	2	1	2	
  chr2	500	1000	ugly	2	+	0	1	1	
  chr3	1000	5000	big	3	-	1	1	0



===========================================================================================
``-both`` Report both the count of hits and the fraction covered from the annotation files
===========================================================================================

.. code-block:: bash

  $ bedtools annotate -both -i variants.bed -files genes.bed conserve.bed known_var.bed
  #chr	start	end	name	score	+/-	cnt1	pct1	cnt2	pct2	cnt3	pct3
  chr1	100	200	nasty	1	-	2	0.500000	1	1.000000	2	0.300000	
  chr2	500	1000	ugly	2	+	0	0.000000	1	0.600000	1	1.000000	
  chr3	1000	5000	big	3	-	1	1.000000	1	0.250000	0	0.000000


  
  
==========================================================================
``-s`` Restrict the reporting to overlaps on the **same** strand.
==========================================================================

.. code-block:: bash

  $ bedtools annotate -s -i variants.bed -files genes.bed conserve.bed known_var.bed
  chr1	100	200	nasty	1	-	0.000000	0.000000	0.000000	
  chr2	500	1000	ugly	2	+	0.000000	0.000000	0.000000	
  chr3	1000	5000	big	3	-	1.000000	0.000000	0.000000



==========================================================================
``-S`` Restrict the reporting to overlaps on the **opposite** strand.
==========================================================================

.. code-block:: bash

  $ bedtools annotate -S -i variants.bed -files genes.bed conserve.bed known_var.bed
  chr1	100	200	nasty	1	-	0.500000	1.000000	0.300000	
  chr2	500	1000	ugly	2	+	0.000000	0.600000	1.000000	
  chr3	1000	5000	big	3	-	0.000000	0.250000	0.000000

