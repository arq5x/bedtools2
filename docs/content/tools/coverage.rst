.. _coverage:

###############
*coverage*
###############
The ``bedtools coverage`` tool computes both the *depth* and *breadth* of coverage of features in file B on the features
in file A. For example, ``bedtools coverage`` can compute the coverage of sequence alignments (file B) across 1
kilobase (arbitrary) windows (file A) tiling a genome of interest. One advantage that ``bedtools coverage``
offers is that it not only *counts* the number of features that overlap an interval in file A, it also
computes the fraction of bases in the interval in A that were overlapped by one or more features. Thus,
``bedtools coverage`` also computes the *breadth* of coverage observed for each interval in A.


.. note::

    If you are trying to compute coverage for very large files and are having trouble
    with excessive memory usage, please presort your data by chromosome and
    then by start position (e.g., ``sort -k1,1 -k2,2n in.bed > in.sorted.bed``
    for BED files) and then use the ``-sorted`` option.  This invokes a 
    memory-efficient algorithm designed for large files.

.. important::

    As of version 2.24.0, the `coverage` tool has changed such that the coverage is
    computed for the A file, not the B file. This changes the command line interface
    to be consistent with the other tools.  Also, the `coverage` tool
    can accept multiple files for the `-b` option. This allows one to measure 
    coverage between a single query (`-a`) file and multiple database files (`-b`) at once!


.. seealso::

    :doc:`../tools/intersect`
    :doc:`../tools/genomecov`
    
===============================
Usage and option summary
===============================
**Usage**:
::

  bedtools coverage [OPTIONS] -a <FILE> \
                               -b <FILE1, FILE2, ..., FILEN>

**(or)**:
::

  coverageBed [OPTIONS] -a <FILE> \
                         -b <FILE1, FILE2, ..., FILEN>


===========================    =========================================================================================================================================================
Option                         Description
===========================    =========================================================================================================================================================
**-a**                         BAM/BED/GFF/VCF file "A". Each feature in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe.
**-b**                         One or more BAM/BED/GFF/VCF file(s) "B". Use "stdin" if passing B with a UNIX pipe.
                               **NEW!!!**: -b may be followed with multiple databases and/or wildcard (*) character(s).
**-abam**                      BAM file A. Each BAM alignment in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe: For example: samtools view -b <BAM> | bedtools intersect -abam stdin -b genes.bed.  **Note**: no longer necessary after version 2.19.0                                                 
**-hist**                      | Report a histogram of coverage for each feature in A as well as a summary histogram for _all_ features in A.
                               | Output (tab delimited) after each feature in A:                 
                               | 1) depth
                               | 2) # bases at depth
                               | 3) size of A
                               | 4) % of A at depth
**-d**                         Report the depth at each position in each A feature. Positions reported are one based. Each position and depth follow the complete A feature.
**-counts**                    Only report the count of overlaps, don't compute fraction, etc. Restricted by -f and -r.
**-f**                         Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
**-F**                         Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).
**-r**                         Require that the fraction of overlap be reciprocal for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
**-e**                         Require that the minimum fraction be satisfied for A _OR_ B. In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of  B is covered. Without -e, both fractions would have to be satisfied.
**-s**                         Force "strandedness". That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
**-S**                         Require different strandedness.  That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.
**-split**                     Treat "split" BAM (i.e., having an "N" CIGAR operation) or BED12 entries as distinct BED intervals.
**-sorted**                    For very large B files, invoke a "sweeping" algorithm that requires position-sorted (e.g., ``sort -k1,1 -k2,2n`` for BED files) input. When using -sorted, memory usage remains low even for very large files.
**-g**                         Specify a genome file the defines the expected chromosome order in the input files for use with the ``-sorted`` option.
**-header**                    Print the header from the A file prior to results.
**-sortout**                   When using *multiple databases* (`-b`), sort the output DB hits for each record.
**-nobuf**                     Disable buffered output. Using this option will cause each line of output to be printed as it is generated, rather than saved in a buffer. This will make printing large output files noticeably slower, but can be useful in conjunction with other software tools and scripts that need to process one line of bedtools output at a time.
**-iobuf**                     Follow with desired integer size of read buffer. Optional suffixes `K/M/G` supported. **Note**: currently has no effect with compressed files.
===========================    =========================================================================================================================================================



==========================================================================
Default behavior
==========================================================================
After each interval in A, ``bedtools coverage`` will report:

1) The number of features in B that overlapped (by at least one base pair) the A interval.
2) The number of bases in A that had non-zero coverage from features in B.
3) The length of the entry in A.
4) The fraction of bases in A that had non-zero coverage from features in B.

Below are the number of features in B (N=...) overlapping A and fraction of bases in A with coverage.

::

  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BED FILE A  ***************     ***************     ******    **************   
  
  BED File B  ^^^^ ^^^^              ^^             ^^^^^^^^^    ^^^ ^^ ^^^^
                ^^^^^^^^                                      ^^^^^ ^^^^^ ^^
  
  Result      [  N=3, 10/15 ]     [  N=1, 2/15 ]     [N=1,6/6]   [N=6, 12/14 ]


For example:

.. code-block:: bash

  $ cat A.bed
  chr1  0   100
  chr1  100 200
  chr2  0   100

  $ cat B.bed
  chr1  10  20
  chr1  20  30
  chr1  30  40
  chr1  100 200

  $ bedtools coverage -a A.bed -b B.bed
  chr1  0   100  3  30  100 0.3000000
  chr1  100 200  1  100 100 1.0000000
  chr2  0   100  0  0   100 0.0000000

  
  
==========================================================================
``-s`` Calculating coverage by strand 
==========================================================================
Use the "**-s**" option if one wants to only count coverage if features in A are on the same strand as the
feature / window in A. This is especially useful for RNA-seq experiments.

For example (note the difference in coverage with and without **-s**:

.. code-block:: bash

  $ cat A.bed
  chr1  0   100 b1  1  +
  chr1  100 200 b2  1  -
  chr2  0   100 b3  1  +

  $ cat B.bed
  chr1  10  20  a1  1  -
  chr1  20  30  a2  1  -
  chr1  30  40  a3  1  -
  chr1  100 200 a4  1  +

  $ bedtools coverage -a A.bed -b B.bed
  chr1  0   100 b1  1  +  3  30  100  0.3000000
  chr1  100 200 b2  1  -  1  100 100  1.0000000
  chr2  0   100 b3  1  +  0  0   100  0.0000000

  $ bedtools coverage -a A.bed -b B.bed -s
  chr1  0   100 b1  1  +  0  0   100  0.0000000
  chr1  100 200 b2  1  -  0  0   100  0.0000000
  chr2  0   100 b3  1  +  0  0   100  0.0000000


==========================================================================
``-hist`` Creating a histogram of coverage for each feature in the A file 
==========================================================================
One should use the "**-hist**" option to create, for each interval in A, a histogram of coverage of the
features in B across A.

In this case, each entire feature in A will be reported, followed by the depth of coverage, the number of
bases at that depth, the size of the feature, and the fraction covered. After all of the features in A have
been reported, a histogram summarizing the coverage among all features in A will be reported.

.. code-block:: bash

  $ cat A.bed
  chr1  0   100 b1  1  +
  chr1  100 200 b2  1  -
  chr2  0   100 b3  1  +

  $ cat B.bed
  chr1  10  20  a1  1  -
  chr1  20  30  a2  1  -
  chr1  30  40  a3  1  -
  chr1  100 200 a4  1  +

  $ bedtools coverage  -a A.bed -b B.bed -hist
  chr1  0   100 b1  1  +  0  70  100  0.7000000
  chr1  0   100 b1  1  +  1  30  100  0.3000000
  chr1  100 200 b2  1  -  1  100 100  1.0000000
  chr2  0   100 b3  1  +  0  100 100  1.0000000
  all   0   170 300 0.5666667
  all   1   130 300 0.4333333


===========================================================================
``-d`` Reporting the per-base of coverage for each feature in the A file 
===========================================================================
One should use the "**-d**" option to create, for each interval in A, a detailed list of coverage at each of the
positions across each A interval.

The output will consist of a line for each one-based position in each A feature, followed by the coverage
detected at that position.

.. code-block:: bash

  $ cat A.bed
  chr1  0  10

  $ cat B.bed
  chr1  0  5
  chr1  3  8
  chr1  4  8
  chr1  5  9

  $ bedtools coverage -a A.bed -b B.bed -d
  chr1  0  10  B  1  1
  chr1  0  10  B  2  1
  chr1  0  10  B  3  1
  chr1  0  10  B  4  2
  chr1  0  10  B  5  3
  chr1  0  10  B  6  3
  chr1  0  10  B  7  3
  chr1  0  10  B  8  3
  chr1  0  10  B  9  1
  chr1  0  10  B  10 0
