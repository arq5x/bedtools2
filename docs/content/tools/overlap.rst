.. _overlap:

###############
*overlap*
###############
**overlap** computes the amount of overlap (in the case of positive values) or distance (in the case of
negative values) between feature coordinates occurring on the same input line and reports the result at
the end of the same line. In this way, it is a useful method for computing custom overlap scores from
the output of other BEDTools.

==========================================================================
Usage and option summary
==========================================================================
Usage:

::

  overlap [OPTIONS] -i <input> -cols s1,e1,s2,e2

===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-i**				             Input file. Use "stdin" for pipes.			 
**-cols**					     Specify the columns (1-based) for the starts and ends of the features for which you'd like to compute the overlap/distance. The columns must be listed in the following order: *start1,end1,start2,end2*
===========================      ===============================================================================================================================================================================================================



==========================================================================
Default behavior
==========================================================================
The default behavior is to compute the amount of overlap between the features you specify based on the
start and end coordinates. For example:

::

  bedtools window -a A.bed -b B.bed -w 10
  chr1  10  20  A  chr1  15  25  B
  chr1  10  20  C  chr1  25  35  D
  
# Now let's say we want to compute the number of base pairs of overlap
# between the overlapping features from the output of windowBed.

::

  bedtools windo -a A.bed -b B.bed -w 10 | overlap -i stdin -cols 2,3,6,7
  chr1  10  20  A  chr1  15  25  B  5
  chr1  10  20  C  chr1  25  35  D  -5

