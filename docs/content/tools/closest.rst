.. _closest:

###############
*closest*
###############


|

.. image:: ../images/tool-glyphs/closest-glyph.png 

|


Similar to :doc:`../tools/intersect`, `closest` searches for overlapping features in A and B. In the event that
no feature in B overlaps the current feature in A, `closest` will report the nearest (that is, least
genomic distance from the start or end of A) feature in B. For example, one might want to find which
is the closest gene to a significant GWAS polymorphism. Note that `closest` will report an
overlapping feature as the closest---that is, it does not restrict to closest *non-overlapping* feature.

.. note::

    ``bedtools closest`` requires that all input files are presorted data by chromosome and
    then by start position (e.g., ``sort -k1,1 -k2,2n in.bed > in.sorted.bed``
    for BED files).

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
**-s**                           | Require same strandedness.  That is, find the closest feature in
                                 | B that overlaps A on the _same_ strand.
                                 | By default, overlaps are reported without respect to strand.

**-S**                           | Require opposite strandedness.  That is, find the closest feature
                                 | in B that overlaps A on the _opposite_ strand.
                                 | By default, overlaps are reported without respect to strand.

**-d**                           | In addition to the closest feature in B,
                                 | report its distance to A as an extra column.
                                 | The reported distance for overlapping features will be 0.

**-D**                           | Like `-d`, report the closest feature in B, and its distance to A
                                 | as an extra column. However unlike `-d`, use negative distances to 
                                 | report upstream features.
                                 | The options for defining which orientation is "upstream" are:
                                 | - `ref`   Report distance with respect to the reference genome.
                                 |           B features with a lower (start, stop) are upstream
                                 | - `a`     Report distance with respect to A.
                                 |           When A is on the - strand, "upstream" means B has a
                                 |           higher (start,stop).
                                 | - `b`     Report distance with respect to B.
                                 |           When B is on the - strand, "upstream" means A has a
                                 |           higher (start,stop).
===========================      ===============================================================================================================================================================================================================




==========================================================================
Default behavior
==========================================================================
**closestBed** first searches for features in B that overlap a feature in A. If overlaps are found, the feature
in B that overlaps the highest fraction of A is reported. If no overlaps are found, **closestBed** looks for
the feature in B that is *closest* (that is, least genomic distance to the start or end of A) to A. For
example, in the figure below, feature B1 would be reported as the closest feature to A1.

::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BED FILE A                             *************
  
  BED File B         ^^^^^^^^                            ^^^^^^
  
  Result                                                 ======
  

For example:

::
  cat A.bed
  chr1  100  200

  cat B.bed
  chr1  500  1000
  chr1  1300 2000

  closestBed -a A.bed -b B.bed
  chr1  100  200  chr1  500  1000



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
