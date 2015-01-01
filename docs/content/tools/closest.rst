.. _closest:

###############
*closest*
###############

Similar to :doc:`../tools/intersect`, `closest` searches for overlapping features in A and B. In the event that
no feature in B overlaps the current feature in A, `closest` will report the nearest (that is, least
genomic distance from the start or end of A) feature in B. For example, one might want to find which
is the closest gene to a significant GWAS polymorphism. Note that `closest` will report an
overlapping feature as the closest---that is, it does not restrict to closest *non-overlapping* feature. The following iconic "cheatsheet" summarizes the funcitonality available through the various optyions provided by the `closest` tool.

|

.. image:: ../images/tool-glyphs/closest-glyph.png 

|



.. note::

    ``bedtools closest`` requires that all input files are presorted data by chromosome and
    then by start position (e.g., ``sort -k1,1 -k2,2n in.bed > in.sorted.bed``
    for BED files).

.. note::

    Reports "none" for chrom and "-1" for all other fields when a feature
    is not found in B on the same chromosome as the feature in A.
    E.g. `none -1  -1`

.. important::

    As of version 2.22.0, the `closest` tool can accept multiple files for
    the `-b` option. This allows one to identify the closest intervals between a single
    query (`-a`) file and multiple database files (`-b`) at once! This functionality
    now requires that all input files be sorted by chromosome and start coordinate
    in an identical manner (e.g., `sort -k1,1 -k2,2n`).


===============================
Usage and option summary
===============================
**Usage**:
::

  bedtools closest [OPTIONS] -a <FILE> \
                             -b <FILE1, FILE2, ..., FILEN>

**(or)**:
::
  
  closestBed [OPTIONS] -a <FILE> \
                       -b <FILE1, FILE2, ..., FILEN>
  

  
===========================      ===============================================================================================================================================================================================================
Option                           Description
===========================      ===============================================================================================================================================================================================================
**-s**                           Require same strandedness.  That is, find the closest feature in B that overlaps A on the _same_ strand. By default, overlaps are reported without respect to strand.

**-S**                           Require opposite strandedness.  That is, find the closest featurein B that overlaps A on the _opposite_ strand. By default, overlaps are reported without respect to strand.

**-d**                           In addition to the closest feature in B, report its distance to A as an extra column. The reported distance for overlapping features will be 0.

**-D**                           | Like `-d`, report the closest feature in B, and its distance to A as an extra column. However unlike `-d`, use negative distances to report upstream features.
                                 | The options for defining which orientation is "upstream" are:
                                 | - `ref`   Report distance with respect to the reference genome.
                                 |           B features with a lower (start, stop) are upstream
                                 | - `a`     Report distance with respect to A.
                                 |           When A is on the - strand, "upstream" means B has a
                                 |           higher (start,stop).
                                 | - `b`     Report distance with respect to B.
                                 |           When B is on the - strand, "upstream" means A has a
                                 |           higher (start,stop).
**-io**                          Ignore features in B that overlap A. That is, we want close, yet not touching features only.

**-iu**                          Ignore features in B that are upstream of features in A. This option requires -D and follows its orientation rules for determining what is "upstream".

**-id**                          Ignore features in B that are downstream of features in A. This option requires -D and follows its orientation rules for determining what is "downstream".

**-t**                           | Specify how ties for closest feature should be handled.  This occurs when two features in B have exactly the same "closeness" with A. By default, all such features in B are reported.
                                 | Here are all the options:
                                 | - `all`    Report all ties (default).
                                 | - `first`  Report the first tie that occurred in the B file.
                                 | - `last`   Report the last tie that occurred in the B file.

**-mdb**                         | Specifiy how multiple databases should be resolved.
                                 | - `each`  Report closest records for each database (default).
                                 | - `all`   Report closest records among all databases.

**-N**                           Require that the query and the closest hit have different names. For BED, the 4th column is compared.

**-header**                      Print the header from the A file prior to results.
===========================      ===============================================================================================================================================================================================================




==========================================================================
Default behavior
==========================================================================
The `closest` tool first searches for features in B that overlap a feature in A. If overlaps are found, the feature in B that overlaps the highest fraction of A is reported. If no overlaps are found, `closestBed` looks for
the feature in B that is *closest* (that is, least genomic distance to the start or end of A) to A. For example, 

For example, consider the case where one of the intervals im B overlaps the interval in B, yet another does not:

.. code-block:: bash

  $ cat a.bed
  chr1  10  20  a1  1 -

  $ cat b.bed
  chr1  7   8   b1  1 -
  chr1  15  25  b2  2 +

  $ bedtools closest -a a.bed -b b.bed
  chr1  10  20  a1  1 - chr1  15  25  b2  2 +


Now compare what happens when neither interval in B overlaps the record in A, yet one is closer than the other.

.. code-block:: bash

  $ cat a.bed
  chr1  10  20  a1  1 -

  $ cat b.bed
  chr1  7   8   b1  1 -
  chr1  30  40  b2  2 +

  $ bedtools closest -a a.bed -b b.bed
  chr1  10  20  a1  1 - chr1  7 8 b1  1

But what if each interval in B is equally close to the interval in A? In this case, the default behavior is to report all intervals in B that are tied for proximity. Check out the `-t` option to adjust this behaviour.

.. code-block:: bash

  $ cat a.bed
  chr1  10  20  a1  1 -

  $ cat b.bed
  chr1  7   8   b1  1 -
  chr1  22  22  b2  2 +

  $ bedtools closest -a a.bed -b b.bed
  chr1  10  20  a1  1 - chr1  7   8   b1  1 -
  chr1  10  20  a1  1 - chr1  22  23  b2  2 +

==========================================================================
``-s`` Enforcing "strandedness" 
==========================================================================
This option behaves the same as the -s option for intersectBed while scanning for the closest
(overlapping or not) feature in B. See the discussion in the intersectBed section for details.



==========================================================================
``-t`` Controlling how ties for "closest" are broken 
==========================================================================
When there are two or more features in B that overlap the *same fraction* of A, **closestBed** will, by
default, report both features in B. Imagine feature A is a SNP and file B contains genes. It can often
occur that two gene annotations (e.g. opposite strands) in B will overlap the SNP. As mentioned, the
default behavior is to report both such genes in B. However, the -t option allows one to optionally
choose the just first or last feature (in terms of where it occurred in the input file, not chromosome
position) that occurred in B.

For example (note the difference between -l 200 and -l 300):

::
  cat A.bed
  chr1  100  101  rs1234

  cat B.bed
  chr1  0  1000  geneA  100  +
  chr1  0  1000  geneB  100  -

  closestBed -a A.bed -b B.bed
  chr1  100  101  rs1234  chr1  0  1000  geneA  100  +
  chr1  100  101  rs1234  chr1  0  1000  geneB  100  -

  closestBed -a A.bed -b B.bed -t all
  chr1  100  101  rs1234  chr1  0  1000  geneA  100  +
  chr1  100  101  rs1234  chr1  0  1000  geneB  100  -

  closestBed -a A.bed -b B.bed -t first
  chr1  100  101  rs1234  chr1  0  1000  geneA  100  +

  closestBed -a A.bed -b B.bed -t last
  chr1  100  101  rs1234  chr1  0  1000  geneB  100  -






==========================================================================
``-d`` Reporting the distance to the closest feature in base pairs 
==========================================================================
ClosestBed will optionally report the distance to the closest feature in the B file using the **-d** option.
When a feature in B overlaps a feature in A, a distance of 0 is reported.

::
  cat A.bed
  chr1  100  200
  chr1  500  600

  cat B.bed
  chr1  500  1000
  chr1  1300 2000

  closestBed -a A.bed -b B.bed -d
  chr1  100  200  chr1  500  1000  300
  chr1  500  600  chr1  500  1000  0
