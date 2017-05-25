.. _sort:

###############
*sort*
###############
The ``bedtools sort`` tool sorts a feature file by chromosome and other criteria.

==========================================================================
Usage and option summary
==========================================================================
**Usage**:
::

  bedtools sort [OPTIONS] -i <BED/GFF/VCF>

**(or)**:
::

  sortBed [OPTIONS] -i <BED/GFF/VCF>

===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-sizeA**                       Sort by feature size in ascending order.
**-sizeD**                       Sort by feature size in descending order.
**-chrThenSizeA**                Sort by chromosome (asc), then by feature size (asc).
**-chrThenSizeD**                Sort by chromosome (asc), then by feature size (desc).
**-chrThenScoreA**               Sort by chromosome (asc), then by score (asc).
**-chrThenScoreD**               Sort by chromosome (asc), then by score (desc).
**-g**                           Define sort order by order of tab-delimited file with chromosome names in the first column.
**-faidx**                       Define sort order by order of tab-delimited file with chromosome names in the first column. Sort by specified chromosome order.
===========================      ===============================================================================================================================================================================================================



==========================================================================
Default behavior
==========================================================================
By default, ``bedtools sort`` sorts a BED file by chromosome and then by start position in ascending order.

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
``bedtools sort`` will also sort a BED file by chromosome and then by other criteria.

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


**Disclaimer:** it should be noted that ``bedtools sort`` is merely a convenience utility, as the UNIX sort utility
will sort BED files more quickly while using less memory. For example, UNIX sort will sort a BED file
by chromosome then by start position in the following manner:

::

  sort -k 1,1 -k2,2n a.bed
  chr1 1   10
  chr1 80  180
  chr1 750 10000
  chr1 800 1000
