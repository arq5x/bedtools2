###############
Example usage
###############
Below are several examples of basic BEDTools usage. Example BED files are provided in the
/data directory of the BEDTools distribution.



==========================================================================
6.1 intersectBed
==========================================================================


6.1.1 Report the base-pair overlap between sequence alignments and genes.
::  
  intersectBed -a reads.bed -b genes.bed



6.1.2 Report whether each alignment overlaps one or more genes. If not, the alignment is not reported.
::
  intersectBed -a reads.bed -b genes.bed -u
  
  

6.1.3 Report those alignments that overlap NO genes. Like "grep -v"
::
  intersectBed -a reads.bed -b genes.bed -v


6.1.4 Report the number of genes that each alignment overlaps.
::
  intersectBed -a reads.bed -b genes.bed -c



6.1.5 Report the entire, original alignment entry for each overlap with a gene.
::
  intersectBed -a reads.bed -b genes.bed -wa



6.1.6 Report the entire, original gene entry for each overlap with a gene.
::
  intersectBed -a reads.bed -b genes.bed -wb
  


6.1.7 Report the entire, original alignment and gene entries for each overlap.
::
  intersectBed -a reads.bed -b genes.bed -wa -wb



6.1.8 Only report an overlap with a repeat if it spans at least 50% of the exon.
::
  intersectBed -a exons.bed -b repeatMasker.bed -f 0.50



6.1.9 Only report an overlap if comprises 50% of the structural variant and 50% of the segmental duplication. Thus, it is reciprocally at least a 50% overlap.
::
  intersectBed -a SV.bed -b segmentalDups.bed -f 0.50 -r

  
  

6.1.10 Read BED A from stdin. For example, find genes that overlap LINEs but not SINEs.
::
  intersectBed -a genes.bed -b LINES.bed | intersectBed -a stdin -b SINEs.bed -v

  
  

6.1.11 Retain only single-end BAM alignments that overlap exons.
::
  intersectBed -abam reads.bam -b exons.bed > reads.touchingExons.bam
  

  
  

6.1.12 Retain only single-end BAM alignments that do not overlap simple sequence
repeats.
::
  intersectBed -abam reads.bam -b SSRs.bed -v > reads.noSSRs.bam



==========================================================================
6.2 pairToBed
==========================================================================



6.2.1 Return all structural variants (in BEDPE format) that overlap with genes on either
end.
::
  pairToBed -a sv.bedpe -b genes > sv.genes



6.2.2 Return all structural variants (in BEDPE format) that overlap with genes on both
end.
::
  pairToBed -a sv.bedpe -b genes -type both > sv.genes


  

6.2.3 Retain only paired-end BAM alignments where neither end overlaps simple
sequence repeats.
::
  pairToBed -abam reads.bam -b SSRs.bed -type neither > reads.noSSRs.bam

  

6.2.4 Retain only paired-end BAM alignments where both ends overlap segmental
duplications.
::
  pairToBed -abam reads.bam -b segdups.bed -type both > reads.SSRs.bam

  
  

6.2.5 Retain only paired-end BAM alignments where neither or one and only one end
overlaps segmental duplications.
::
  pairToBed -abam reads.bam -b segdups.bed -type notboth > reads.notbothSSRs.bam


  
  
  
  
==========================================================================
6.3 pairToPair
==========================================================================


6.3.1 Find all SVs (in BEDPE format) in sample 1 that are also in sample 2.
::
  pairToPair -a 1.sv.bedpe -b 2.sv.bedpe | cut -f 1-10 > 1.sv.in2.bedpe



6.3.2 Find all SVs (in BEDPE format) in sample 1 that are not in sample 2.
::
  pairToPair -a 1.sv.bedpe -b 2.sv.bedpe -type neither | cut -f 1-10 >
1.sv.notin2.bedpe





==========================================================================
6.4 bamToBed
==========================================================================


6.4.1 Convert BAM alignments to BED format.
::
  bamToBed -i reads.bam > reads.bed


6.4.2 Convert BAM alignments to BED format using the BAM edit distance (NM) as the
BED "score".
::
  bamToBed -i reads.bam -ed > reads.bed


6.4.3 Convert BAM alignments to BEDPE format.
::
  bamToBed -i reads.bam -bedpe > reads.bedpe
  
  

  

==========================================================================
6.5 windowBed
==========================================================================



6.5.1 Report all genes that are within 10000 bp upstream or downstream of CNVs.
::
  windowBed -a CNVs.bed -b genes.bed -w 10000



6.5.2 Report all genes that are within 10000 bp upstream or 5000 bp downstream of
CNVs.
::
  windowBed -a CNVs.bed -b genes.bed -l 10000 -r 5000


6.5.3 Report all SNPs that are within 5000 bp upstream or 1000 bp downstream of genes.
Define upstream and downstream based on strand.
::
  windowBed -a genes.bed -b snps.bed -l 5000 -r 1000 -sw

  
  
  
  
==========================================================================
6.6 closestBed
==========================================================================
Note: By default, if there is a tie for closest, all ties will be reported. **closestBed** allows overlapping
features to be the closest.



6.6.1 Find the closest ALU to each gene.
::
  closestBed -a genes.bed -b ALUs.bed


6.6.2 Find the closest ALU to each gene, choosing the first ALU in the file if there is a
tie.
::
  closestBed -a genes.bed -b ALUs.bed -t first



6.6.3 Find the closest ALU to each gene, choosing the last ALU in the file if there is a
tie.
::
  closestBed -a genes.bed -b ALUs.bed -t last


  
  
  
==========================================================================
6.7 subtractBed
==========================================================================
Note: If a feature in A is entirely "spanned" by any feature in B, it will not be reported.



6.7.1 Remove introns from gene features. Exons will (should) be reported.
::
  subtractBed -a genes.bed -b introns.bed
  
  
==========================================================================
6.8 mergeBed
==========================================================================


6.8.1 Merge overlapping repetitive elements into a single entry.
::
  mergeBed -i repeatMasker.bed



6.8.2 Merge overlapping repetitive elements into a single entry, returning the number of
entries merged.
::
  mergeBed -i repeatMasker.bed -n


6.8.3 Merge nearby (within 1000 bp) repetitive elements into a single entry.
::
  mergeBed -i repeatMasker.bed -d 1000


==========================================================================
6.9 coverageBed
==========================================================================


6.9.1 Compute the coverage of aligned sequences on 10 kilobase "windows" spanning the
genome.
::
  coverageBed -a reads.bed -b windows10kb.bed | head
  chr1 0     10000 0  10000 0.00
  chr1 10001 20000 33 10000 0.21
  chr1 20001 30000 42 10000 0.29
  chr1 30001 40000 71 10000 0.36

  

6.9.2 Compute the coverage of aligned sequences on 10 kilobase "windows" spanning the
genome and created a BEDGRAPH of the number of aligned reads in each window for
display on the UCSC browser.
::
  coverageBed -a reads.bed -b windows10kb.bed | cut -f 1-4 > windows10kb.cov.bedg

  

6.9.3 Compute the coverage of aligned sequences on 10 kilobase "windows" spanning the
genome and created a BEDGRAPH of the fraction of each window covered by at least
one aligned read for display on the UCSC browser.
::
  coverageBed -a reads.bed -b windows10kb.bed | awk ¡®{OFS="\t"; print $1,$2,$3,$6}¡¯
  > windows10kb.pctcov.bedg




==========================================================================
6.10 complementBed
==========================================================================


6.10.1 Report all intervals in the human genome that are not covered by repetitive
elements.
::
  complementBed -i repeatMasker.bed -g hg18.genome


  
==========================================================================
6.11 shuffleBed
==========================================================================


6.11.1 Randomly place all discovered variants in the genome. However, prevent them
from being placed in know genome gaps.
::
  shuffleBed -i variants.bed -g hg18.genome -excl genome_gaps.bed


6.11.2 Randomly place all discovered variants in the genome. However, prevent them
from being placed in know genome gaps and require that the variants be randomly
placed on the same chromosome.
::
  shuffleBed -i variants.bed -g hg18.genome -excl genome_gaps.bed -chrom
