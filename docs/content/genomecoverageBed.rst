###############
5.10 genomeCoverageBed
###############
**genomeCoverageBed** computes a histogram of feature coverage (e.g., aligned sequences) for a given
genome. Optionally, by using the **-d** option, it will report the depth of coverage at *each base* on each
chromosome in the genome file (**-g**).

==========================================================================
5.10.1 Usage and option summary
==========================================================================
Usage:
::
  genomeCoverageBed [OPTIONS] -i <BED> -g <GENOME>
  
NOTE: genomeCoverageBed requires that the input BED file be sorted by
chromosome. A simple sort -k1,1 will suffice.

===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-ibam**				         BAM file as input for coverage. Each BAM alignment in A added to the total coverage for the genome. Use "stdin" if passing it with a UNIX pipe: For example:
                                 | samtools view -b <BAM> | genomeCoverageBed -ibam stdin -g hg18.genome								 
**-d**					         Report the depth at each genome position. *Default behavior is to report a histogram*.
**-max**                         Combine all positions with a depth >= max into a single bin in the histogram.
**-bg**                          Report depth in BedGraph format. For details, see: http://genome.ucsc.edu/goldenPath/help/bedgraph.html
**-bga**                         Report depth in BedGraph format, as above (i.e., -bg). However with this option, regions with zero coverage are also reported. This allows one to quickly extract all regions of a genome with 0 coverage by applying: "grep -w 0$" to the output.
**-split**                       Treat "split" BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR "N" and "D" operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12).
**-strand**                      Calculate coverage of intervals from a specific strand. With BED files, requires at least 6 columns (strand is column 6).
===========================      ===============================================================================================================================================================================================================




==========================================================================
5.10.2 Default behavior
==========================================================================
By default, **genomeCoverageBed** will compute a histogram of coverage for the genome file provided.
The default output format is as follows:
1. chromosome (or entire genome)
2. depth of coverage from features in input file
3. number of bases on chromosome (or genome) with depth equal to column 2.
4. size of chromosome (or entire genome) in base pairs
5. fraction of bases on chromosome (or entire genome) with depth equal to column 2.

For example:
::
  cat A.bed
  chr1  10  20
  chr1  20  30
  chr2  0   500

  cat my.genome
  chr1  1000
  chr2  500

  genomeCoverageBed -i A.bed -g my.genome
  chr1   0  980  1000  0.98
  chr1   1  20   1000  0.02
  chr2   1  500  500   1
  genome 0  980  1500  0.653333
  genome 1  520  1500  0.346667

  
  
  
==========================================================================
5.10.3 (-max)Controlling the histogram's maximum depth 
==========================================================================
Using the **-max** option, **genomeCoverageBed** will "lump" all positions in the genome having feature
coverage greather than or equal to **max** into the **max** histogram bin. For example, if one sets **-max**
equal to 50, the max depth reported in the output will be 50 and all positions with a depth >= 50 will
be represented in bin 50.

==========================================================================
5.10.4 (-d)Reporting "per-base" genome coverage 
==========================================================================
Using the **-d** option, **genomeCoverageBed** will compute the depth of feature coverage for each base
on each chromosome in genome file provided.

The "per-base" output format is as follows:
1. chromosome
2. chromosome position
3. depth (number) of features overlapping this chromosome position.

For example:
::
  cat A.bed
  chr1  10  20
  chr1  20  30
  chr2  0   500

  cat my.genome
  chr1  1000
  chr2  500

  genomeCoverageBed -i A.bed -g my.genome -d | head -15 | tail -n 10
  chr1  6  0
  chr1  7  0
  chr1  8  0
  chr1  9  0
  chr1  10 0
  chr1  11 1
  chr1  12 1
  chr1  13 1
  chr1  14 1
  chr1  15 1

  
  
==========================================================================
5.1.13 (-split)Reporting coverage with spliced alignments or blocked BED features 
==========================================================================
As described in section 1.3.19, genomeCoverageBed will, by default, screen for overlaps against the
entire span of a spliced/split BAM alignment or blocked BED12 feature. When dealing with RNA-seq
reads, for example, one typically wants to only screen for overlaps for the portions of the reads that
come from exons (and ignore the interstitial intron sequence). The **-split** command allows for such
overlaps to be performed.

For additional details, please visit the Usage From The Wild site and have a look at example 5,
contributed by Assaf Gordon.


