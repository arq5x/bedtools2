.. _multicov:

###############
*multicov*
###############
``bedtools multicov``, reports the count of alignments from multiple 
position-sorted and indexed BAM files that overlap intervals in a BED file.
Specifically, for each BED interval provided, it reports a separate count of
overlapping alignments from each BAM file.

.. note::

    ``bedtools multicov`` depends upon index BAM files in order to count the
    number of overlaps in each BAM file.  As such, each BAM file should be
    position sorted (``samtool sort aln.bam aln.sort``) and 
    indexed (``samtools index aln.sort.bam``) with either samtools or bamtools.

    
==========================================================================
Usage and option summary
==========================================================================
**Usage**:
::

  bedtools multicov [OPTIONS] -bams BAM1 BAM2 BAM3 ... BAMn -bed  <BED/GFF/VCF>

**(or)**:
::

  multiBamCov [OPTIONS] -bams BAM1 BAM2 BAM3 ... BAMn -bed  <BED/GFF/VCF>
  
  
============    ================================================================
 Option          Description
============    ================================================================
**-split**      Treat "split" BAM or BED12 entries as distinct BED intervals.
**-s**          Require same strandedness.  That is, only report hits in B
                that overlap A on the _same_ strand. By default, overlaps are 
                reported without respect to strand.
**-S**          Require different strandedness. That is, only report hits in B
                that overlap A on the _opposite_ strand. By default, overlaps 
                are reported without respect to strand.
**-f**          Minimum overlap required as a fraction of each -bed record. Default is 
                1E-9 (i.e., 1bp).
**-r**          Require that the fraction overlap is reciprocal for -bed and B. In 
                other words, if -f is 0.90 and -r is used, this requires that 
                B overlap 90% of the -bed record and the -bed record _also_ overlaps 90% of B.
**-q**          Minimum mapping quality (MAPQ) allowed. Default is 0.
**-D**          Include duplicate reads.  Default counts non-duplicates only
**-F**          Include failed-QC reads.  Default counts pass-QC reads only
**-p**          Only count proper pairs.  Default counts all alignments with
                ``MAPQ > -q`` argument, regardless of the BAM FLAG field.
============    ================================================================


==========================================================================
Default behavior.
==========================================================================
By default, ``multicov`` will report the count of alignments in each input
BAM file that overlap.

.. code-block:: bash

   $ cat ivls-of-interest.bed
   chr1 0   10000   ivl1
   chr1 10000   20000   ivl2
   chr1 20000   30000   ivl3
   chr1 30000   40000   ivl4
   
   $ bedtools multicov -bams aln1.bam aln2.bam aln3.bam -bed ivls-of-interest.bed
   chr1	0	10000	ivl1	100 2234    0
   chr1	10000	20000	ivl2	123 3245    1000
   chr1	20000	30000	ivl3	213 2332    2034
   chr1	30000	40000	ivl4	335 7654    0


The output of ``multicov`` reflects a distinct report of the overlapping
alignments for each record in the ``-bed`` file.  In the example above, each 
line of the output reflects **a)** the original line from the ``-bed`` file 
followed by **b)** the count of alignments that overlap the ``-bed`` interval
from each input ``-bam`` file.  In the example above, the output consists of
7 columns: the first four of which are the columns from the ``-bed`` file and
the last 3 are the count of overlapping alignments from the 3 input ``-bam`` 
files.  The order of the counts reflects the order of the files given on the 
command line.

.. note::

    ``bedtools multicov`` will work with a single BAM as well.

.. code-block:: bash

   $ bedtools multicov -bams aln1.bam -bed ivls-of-interest.bed
   chr1	0	10000	ivl1	100
   chr1	10000	20000	ivl2	123
   chr1	20000	30000	ivl3	213
   chr1	30000	40000	ivl4	335
   
