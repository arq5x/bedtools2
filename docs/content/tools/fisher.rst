.. _fisher:

########
*fisher*
########

Perform fisher's exact test on the non/overlap between 2 files.

|

Traditionally, in order to test whether 2 sets of intervals are related
spatially, we resort to shuffling the genome and checking the simulated
(shuffled) versus the observed. We can do the same analytically for many
scenarios using 
`Fisher's Exact Test`_ .



Given a pair of input files `-a` and `-b` in the usual BedTools parlance:

.. code-block:: bash

  $ cat a.bed
  chr1  10  20
  chr1  30  40

  $ cat b.bed
  chr1  15   25

And a genome of 500 bases:

.. code-block:: bash

    $ echo -e "chr1\t500" > t.genome

We may wish to know **if the amount of overlap between the 2 sets of intervals is
more than we would expect given their coverage and the size of the genome**. We
can do this with ``fisher`` as:

.. code-block:: bash

    $ bedtools fisher -a a.bed -b b.bed -g t.genome
    # Contingency Table
    #_________________________________________
    #           | not in -b    | in -b        |
    # not in -a | 475          | 5           |
    #     in -a | 15           | 5            |
    #_________________________________________
    # p-values for fisher's exact test
    left    right   two-tail    ratio
    1       1.3466e-05      1.3466e-05      31.667


Where we can see the constructed contingency table and the pvalues for left, right
and two-tail tests.
From here, we can say that given **500 bases** of genome, it is unlikely that a region of
20 bases total from `-a` and 10 bases total from `-b` would share 5 bases if the regions
were randomly distributed. *Consult your statistician for a more precise definition*.

The above was highly significant, but if our genome were only **50 bases**:

.. code-block:: bash

    $ echo -e "chr1\t50" > t.genome
    $ bedtools fisher -a a.bed -b b.bed -g t.genome
    # Contingency Table
    #_________________________________________
    #           | not in -b    | in -b        |
    # not in -a | 25           | 15           |
    #     in -a | 5            | 5            |
    #_________________________________________
    # p-values for fisher's exact test
    left    right   two-tail    ratio
    0.86011 0.35497 0.49401 1.667


We can see that neither tail is significant. Intuitively, this makes sense; 
if we randomly place 20 (from `-a`), and 10 (from `-b`) bases of intervals
within 50 bases, it doesn't seem unlikely that we'd see 5 bases of overlap.


.. note::

    The ``fisher`` tool requires that your data is pre-sorted by chromosome and
    then by start position (e.g., ``sort -k1,1 -k2,2n in.bed > in.sorted.bed``
    for BED files).

    This uses Heng Li's implementation of Fisher's exact test in kfunc.c.

.. seealso::

    :doc:`../tools/jaccard`
    :doc:`../tools/reldist`
    :doc:`../tools/intersect`
    

===============================
Usage and option summary
===============================
**Usage**:
::

  bedtools fisher [OPTIONS] -a <BED/GFF/VCF> -b <BED/GFF/VCF> -g <genome>


===========================    =========================================================================================================================================================
Option                         Description
===========================    =========================================================================================================================================================
**-a**                           BED/GFF/VCF file A. Each feature in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe.
**-b**                           BED/GFF/VCF file B. Use "stdin" if passing B with a UNIX pipe.
**-g**                           genome file listing chromosome size.
**-f**                           Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
**-r**                           Require that the fraction of overlap be reciprocal for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
**-s**                         Force "strandedness". That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
**-S**                         Require different strandedness.  That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.
**-split**                     Treat "split" BAM (i.e., having an "N" CIGAR operation) or BED12 entries as distinct BED intervals.
===========================    =========================================================================================================================================================

.. _Fisher's Exact Test: http://en.wikipedia.org/wiki/Fisher's_exact_test
