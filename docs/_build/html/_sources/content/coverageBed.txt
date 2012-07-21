###############
5.9 coverageBed
###############
**coverageBed** computes both the *depth* and *breadth* of coverage of features in file A across the features
in file B. For example, **coverageBed** can compute the coverage of sequence alignments (file A) across 1
kilobase (arbitrary) windows (file B) tiling a genome of interest. One advantage that **coverageBed**
offers is that it not only *counts* the number of features that overlap an interval in file B, it also
computes the fraction of bases in B interval that were overlapped by one or more features. Thus,
**coverageBed** also computes the *breadth* of coverage for each interval in B.

==========================================================================
5.9.1 Usage and option summary
==========================================================================
Usage:
::
  coverageBed [OPTIONS] -a <BED/GFF/VCF> -b <BED/GFF/VCF>
  
===========================      ===============================================================================================================================================================================================================
Option                           Description
===========================      ===============================================================================================================================================================================================================
**-abam**				         BAM file A. Each BAM alignment in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe: For example:

                                 | samtools view -b <BAM> | intersectBed -abam stdin -b genes.bed
								 
**-s**					         Force strandedness. That is, only features in A are only counted towards coverage in B if they are the same strand. *By default, this is disabled and coverage is counted without respect to strand*.
**-hist**                        Report a histogram of coverage for each feature in B as well as a summary histogram for _all_ features in B.

                                 | Output (tab delimited) after each feature in B:
								 
								 | 1) depth
								 | 2) # bases at depth
								 | 3) size of B
								 | 4) % of B at depth
**-d**                           Report the depth at each position in each B feature. Positions reported are one based. Each position and depth follow the complete B feature.
**-split**                       Treat "split" BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR "N" and "D" operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12).
===========================      ===============================================================================================================================================================================================================






==========================================================================
5.9.2 Default behavior
==========================================================================
After each interval in B, **coverageBed** will report:

1) The number of features in A that overlapped (by at least one base pair) the B interval.
2) The number of bases in B that had non-zero coverage from features in A.
3) The length of the entry in B.
4) The fraction of bases in B that had non-zero coverage from features in A.

Below are the number of features in A (N=...) overlapping B and fraction of bases in B with coverage.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BED FILE B  ***************     ***************     ******    **************   
  
  BED File A  ^^^^ ^^^^              ^^             ^^^^^^^^^    ^^^ ^^ ^^^^
                ^^^^^^^^                                      ^^^^^ ^^^^^ ^^
  
  Result      [  N=3, 10/15 ]     [  N=1, 2/16 ]     [N=1,6/6]   [N=5, 11/12 ]


For example:
::
  cat A.bed
  chr1  10  20
  chr1  20  30
  chr1  30  40
  chr1  100 200

  cat B.bed
  chr1  0   100
  chr1  100 200
  chr2  0   100

  coverageBed -a A.bed -b B.bed
  chr1  0   100  3  30  100 0.3000000
  chr1  100 200  1  100 100 1.0000000
  chr2  0   100  0  0   100 0.0000000

  
  
==========================================================================
5.9.4 (-s)Calculating coverage by strand 
==========================================================================
Use the "**-s**" option if one wants to only count coverage if features in A are on the same strand as the
feature / window in B. This is especially useful for RNA-seq experiments.

For example (note the difference in coverage with and without **-s**:
::
  cat A.bed
  chr1  10  20  a1  1  -
  chr1  20  30  a2  1  -
  chr1  30  40  a3  1  -
  chr1  100 200 a4  1  +

  cat B.bed
  chr1  0   100 b1  1  +
  chr1  100 200 b2  1  -
  chr2  0   100 b3  1  +

  coverageBed -a A.bed -b B.bed
  chr1  0   100 b1  1  +  3  30  100  0.3000000
  chr1  100 200 b2  1  -  1  100 100  1.0000000
  chr2  0   100 b3  1  +  0  0   100  0.0000000

  coverageBed -a A.bed -b B.bed -s
  chr1  0   100 b1  1  +  0  0   100  0.0000000
  chr1  100 200 b2  1  -  0  0   100  0.0000000
  chr2  0   100 b3  1  +  0  0   100  0.0000000

==========================================================================
5.9.5 (-hist)Creating a histogram of coverage for each feature in the B file 
==========================================================================
One should use the "**-hist**" option to create, for each interval in B, a histogram of coverage of the
features in A across B.

In this case, each entire feature in B will be reported, followed by the depth of coverage, the number of
bases at that depth, the size of the feature, and the fraction covered. After all of the features in B have
been reported, a histogram summarizing the coverage among all features in B will be reported.
::
  cat A.bed
  chr1  10  20  a1  1  -
  chr1  20  30  a2  1  -
  chr1  30  40  a3  1  -
  chr1  100 200 a4  1  +

  cat B.bed
  chr1  0   100 b1  1  +
  chr1  100 200 b2  1  -
  chr2  0   100 b3  1  +

  coverageBed -a A.bed -b B.bed -hist
  chr1  0   100 b1  1  +  0  70  100  0.7000000
  chr1  0   100 b1  1  +  1  30  100  0.3000000
  chr1  100 200 b2  1  -  1  100 100  1.0000000
  chr2  0   100 b3  1  +  0  100 100  1.0000000
  all   0   170 300 0.5666667
  all   1   130 300 0.4333333



==========================================================================
5.9.6 (-hist)Reporting the per-base of coverage for each feature in the B file 
==========================================================================
One should use the "**-d**" option to create, for each interval in B, a detailed list of coverage at each of the
positions across each B interval.

The output will consist of a line for each one-based position in each B feature, followed by the coverage
detected at that position.
::
  cat A.bed
  chr1  0  5
  chr1  3  8
  chr1  4  8
  chr1  5  9

  cat B.bed
  chr1  0  10

  coverageBed -a A.bed -b B.bed -d
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

  
  
==========================================================================
5.9.7 (-split)Reporting coverage with spliced alignments or blocked BED features 
==========================================================================
As described in section 1.3.19, coverageBed will, by default, screen for overlaps against the entire span
of a spliced/split BAM alignment or blocked BED12 feature. When dealing with RNA-seq reads, for
example, one typically wants to only tabulate coverage for the portions of the reads that come from
exons (and ignore the interstitial intron sequence). The **-split** command allows for such coverage to be
performed.
