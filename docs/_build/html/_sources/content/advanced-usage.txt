###############
Advanced usage
###############


==========================================================================
7.1 Mask all regions in a genome except for targeted capture regions.
==========================================================================
# Add 500 bp up and downstream of each probe
::
  slopBed -i probes.bed -b 500 > probes.500bp.bed
  
# Get a BED file of all regions not covered by the probes (+500 bp up/down)
::
  complementBed -i probes.500bp.bed -g hg18.genome > probes.500bp.complement.bed
  
# Create a masked genome where all bases are masked except for the probes +500bp
::
  maskFastaFromBed -in hg18.fa -bed probes.500bp.complement.bed -fo hg18.probecomplement.
  masked.fa


==========================================================================
7.2 Screening for novel SNPs.
==========================================================================
# Find all SNPs that are not in dbSnp and not in the latest 1000 genomes calls
::
  intersectBed -a snp.calls.bed -b dbSnp.bed -v | intersectBed -a stdin -b 1KG.bed
  -v > snp.calls.novel.bed



==========================================================================
7.3 Computing the coverage of features that align entirely within an
interval.
==========================================================================
# By default, coverageBed counts any feature in A that overlaps B by >= 1 bp. If
you want to require that a feature align entirely within B for it to be counted,
you can first use intersectBed with the "-f 1.0" option.
::
  intersectBed -a features.bed -b windows.bed -f 1.0 | coverageBed -a stdin -b
  windows.bed > windows.bed.coverage


==========================================================================
7.4 Computing the coverage of BAM alignments on exons.
==========================================================================
# One can combine SAMtools with BEDtools to compute coverage directly from the BAM
data by using bamToBed.
::
  bamToBed -i reads.bam | coverageBed -a stdin -b exons.bed > exons.bed.coverage
  
# Take it a step further and require that coverage be from properly-paired reads.
::
  samtools view -bf 0x2 reads.bam | bamToBed -i stdin | coverageBed -a stdin -b
  exons.bed > exons.bed.proper.coverage



==========================================================================
7.5 Computing coverage separately for each strand.
==========================================================================
# Use grep to only look at forward strand features (i.e. those that end in "+").
::
  bamToBed -i reads.bam | grep \+$ | coverageBed -a stdin -b genes.bed >
  genes.bed.forward.coverage

# Use grep to only look at reverse strand features (i.e. those that end in "-").
::
  bamToBed -i reads.bam | grep \-$ | coverageBed -a stdin -b genes.bed >
  genes.bed.forward.coverage


  
==========================================================================
7.6 Find structural variant calls that are private to one sample.
==========================================================================
# :
::
  pairToPair -a sample1.sv.bedpe -b othersamples.sv.bedpe -type neither >
  sample1.sv.private.bedpe
  
  

==================================================================================
7.7 Exclude SV deletions that appear to be ALU insertions in the reference genome.
==================================================================================
# We'll require that 90% of the inner span of the deletion be overlapped by a
recent ALU.
::
  pairToBed -a deletions.sv.bedpe -b ALUs.recent.bed -type notispan -f 0.80 >
  deletions.notALUsinRef.bedpe