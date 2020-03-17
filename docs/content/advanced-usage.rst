###############
Advanced usage
###############


==========================================================================
Mask all regions in a genome except for targeted capture regions.
==========================================================================

Step 1. Add 500 bp up and downstream of each probe

.. code-block:: bash

  bedtools slop -i probes.bed -g hg18.genome -b 500 > probes.500bp.bed
  
NB genome is two column chromosome size list - i.e. https://genome.ucsc.edu/goldenpath/help/hg18.chrom.sizes

Step 2. Get a BED file of all regions not covered by the probes (+500 bp up/down)


.. code-block:: bash

  bedtools complement -i probes.500bp.bed -g hg18.genome > probes.500bp.complement.bed
  

Step 3. Create a masked genome where all bases are masked except for the probes +500bp

.. code-block:: bash

  bedtools maskfasta -fi hg18.fa -bed probes.500bp.complement.bed \
  -fo hg18.probecomplement.masked.fa



==========================================================================
Screening for novel SNPs.
==========================================================================
Find all SNPs that are not in dbSnp and not in the latest 1000 genomes calls

.. code-block:: bash

  bedtools intersect -a snp.calls.bed -b dbSnp.bed -v | \ 
  bedtools intersect -a - -b 1KG.bed -v | \
  > snp.calls.novel.bed


==========================================================================
Computing the coverage of features that align entirely within an interval.
==========================================================================

By default, bedtools ``coverage`` counts any feature in A that overlaps B 
by >= 1 bp. If you want to require that a feature align entirely within B for 
it to be counted, you can first use intersectBed with the "-f 1.0" option.

.. code-block:: bash

  bedtools intersect -a features.bed -b windows.bed -f 1.0 | \
  bedtools coverage -a windows.bed -b - \
  > windows.bed.coverage


==========================================================================
Computing the coverage of BAM alignments on exons.
==========================================================================
One can combine ``samtools`` with ``bedtools`` to compute coverage directly 
from the BAM data by using ``bamtobed``.

.. code-block:: bash

  bedtools bamtobed -i reads.bam | \
  bedtools coverage -a exons.bed -b - \
  > exons.bed.coverage
  

Take it a step further and require that coverage be from properly-paired reads.

.. code-block:: bash

  samtools view -uf 0x2 reads.bam | \
  coverageBed -abam - -b exons.bed \
  > exons.bed.proper.coverage



==========================================================================
Computing coverage separately for each strand.
==========================================================================
Use grep to only look at forward strand features (i.e. those that end in "+").

.. code-block:: bash

  bedtools bamtobed -i reads.bam | \
  grep \+$  | \
  bedtools coverage -a - -b genes.bed \
  > genes.bed.forward.coverage

Use grep to only look at reverse strand features (i.e. those that end in "-").

.. code-block:: bash

  bedtools bamtobed -i reads.bam | \
  grep \-$ | \
  bedtools coverage -a - -b genes.bed \
  > genes.bed.reverse.coverage

