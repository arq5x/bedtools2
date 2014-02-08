.. _sort:

###############
*sort*
###############
**sortBed** sorts a feature file by chromosome and other criteria.

==========================================================================
Usage and option summary
==========================================================================
Usage:

::
  sortBed [OPTIONS] -i <BED/GFF/VCF>

===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-sizeA**				         Sort by feature size in ascending order.					 
**-sizeD**					     Sort by feature size in descending order.
**-chrThenSizeA**                Sort by chromosome, then by feature size (asc).
**-chrThenSizeD**                Sort by chromosome, then by feature size (desc).
**-chrThenScoreA**               Sort by chromosome, then by score (asc).
**-chrThenScoreD**               Sort by chromosome, then by score (desc).
===========================      ===============================================================================================================================================================================================================



==========================================================================
Default behavior
==========================================================================
By default, **sortBed** sorts a BED file by chromosome and then by start position in ascending order.

For example:

::
  cat A.bed
  chr1 800 1000
  chr1 80  180
  chr1 1   10
  chr1 750 10000

  sortBed -i A.bed
  chr1 1   10
  chr1 80  180
  chr1 750 10000
  chr1 800 1000


  
  
==========================================================================
Optional sorting behavior
==========================================================================
**sortBed** will also sorts a BED file by chromosome and then by other criteria.

For example, to sort by chromosome and then by feature size (in descending order):

::
  cat A.bed
  chr1 800 1000
  chr1 80  180
  chr1 1   10
  chr1 750 10000

  sortBed -i A.bed -sizeD
  chr1 750 10000
  chr1 800 1000
  chr1 80  180
  chr1 1   10
  

**Disclaimer:** it should be noted that **sortBed** is merely a convenience utility, as the UNIX sort utility
will sort BED files more quickly while using less memory. For example, UNIX sort will sort a BED file
by chromosome then by start position in the following manner:

::
  sort -k 1,1 -k2,2n a.bed
  chr1 1   10
  chr1 80  180
  chr1 750 10000
  chr1 800 1000

