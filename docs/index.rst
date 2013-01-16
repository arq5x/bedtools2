================================================================
**bedtools**: *a powerful toolset for genome arithmetic*
================================================================

Collectively, the **bedtools** utilities are a swiss-army knife of tools
for a wide-range of genomics analysis tasks. The most widely-used
tools enable *genome arithmetic*: that is, set theory on the genome.  For 
example, **bedtools** allows one to *intersect*, *merge*, *count*, *complement*,
and *shuffle* genomic intervals from multiple files in widely-used 
genomic file formats such as BAM, BED, GFF/GTF, VCF. 

While each individual tool is designed to do a relatively simple task (e.g., 
*intersect* two interval files), quite sophisticated analyses can be conducted
by combining multiple bedtools operations on the UNIX command line.

=================
Table of contents
=================
.. toctree::
   :maxdepth: 1
   :numbered:

   content/overview
   content/installation
   content/quick-start
   content/general-usage
   content/bedtools-suite
   content/example-usage
   content/advanced-usage
   content/tips-and-tricks
   content/faq
   content/related-tools
   

=================
Brief example
=================
Let's imagine you have a BED file of ChiP-seq peaks from two different
experiments. You want to identify peaks that were observed in *both* experiments
(requiring 50% reciprocal overlap) and for those peaks, you want to find to 
find the closest, non-overlapping gene. Such an analysis could be conducted 
with two, relatively simple bedtools commands.

.. code-block:: bash

    # intersect the peaks from both experiments.
    # -f 0.50 combined with -r requires 50% reciprocal overlap between the 
    # peaks from each experiment.
    $ bedtools intersect -a exp1.bed -b exp2.bed -f 0.50 -r > both.bed
    
    # find the closest, non-overlapping gene for each interval where
    # both experiments had a peak
    # -io ignores overlapping intervals and returns only the closest, 
    # non-overlapping interval (in this case, genes)
    $ bedtools closest -a both.bed -b genes.bed -io > both.nearest.genes.txt

==========
License
==========
bedtools is freely available under a GNU Public License (Version 2).

=====================================
Acknowledgments
=====================================

To do.
    

=================
Mailing list
=================
If you have questions, requests, or bugs to report, please email the
`bedtools mailing list <https://groups.google.com/forum/?fromgroups#!forum/bedtools-discuss>`_

