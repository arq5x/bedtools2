#########################################
5.1 intersect
#########################################

By far, the most common question asked of two sets of genomic features is whether or not any of the
features in the two sets "overlap" with one another. This is known as feature intersection. **bedtools intersect**
allows one to screen for overlaps between two sets of genomic features. Moreover, it allows one to have
fine control as to how the intersections are reported. **bedtools intersect** works with both BED/GFF/VCF
and BAM files as input.

===============================
5.1.1 Usage and option summary
===============================
**Usage**:
::

  bedtools intersect [OPTIONS] [-a <BED/GFF/VCF> || -abam <BAM>] -b <BED/GFF/VCF>
  
  intersectBed [OPTIONS] [-a <BED/GFF/VCF> || -abam <BAM>] -b <BED/GFF/VCF>
  
  

===========================      =========================================================================================================================================================
Option                           Description
===========================      =========================================================================================================================================================
**-a**				             BED/GFF/VCF file A. Each feature in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe.
**-b**					         BED/GFF/VCF file B. Use "stdin" if passing B with a UNIX pipe.
**-abam**					     BAM file A. Each BAM alignment in A is compared to B in search of overlaps. Use "stdin" if passing A with a UNIX pipe: For example: samtools view -b <BAM> | bedtools intersect -abam stdin -b genes.bed                                                   
**-ubam**					     Write uncompressed BAM output. The default is write compressed BAM output.
**-bed**					     When using BAM input (-abam), write output as BED. The default is to write output in BAM when using -abam. For example:   bedtools intersect -abam reads.bam -b genes.bed -bed                              
**-wa**					         Write the original entry in A for each overlap.
**-wb** 				         Write the original entry in B for each overlap. Useful for knowing what A overlaps. Restricted by -f and -r.
**-wo** 				         Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r.
**-wao** 						 Write the original A and B entries plus the number of base pairs of overlap between the two features. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0. Restricted by -f and -r.
**-u**						     Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B. Restricted by -f and -r.
**-c** 			                 For each entry in A, report the number of hits in B while restricting to -f. Reports 0 for A entries that have no overlap with B. Restricted by -f and -r.
**-v**	 			             Only report those entries in A that have no overlap in B. Restricted by -f and -r.
**-f**					         Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
**-r**						     Require that the fraction of overlap be reciprocal for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
**-s**						     Force "strandedness". That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
**-split**					     Treat "split" BAM (i.e., having an "N" CIGAR operation) or BED12 entries as distinct BED intervals.
===========================      =========================================================================================================================================================


===============================
5.1.2 Default behavior
===============================
By default, if an overlap is found, **bedtools intersect** reports the shared interval between the two
overlapping features.

For example:
::
  cat A.bed
  chr1  10  20
  chr1  30  40

  cat B.bed
  chr1  15   20

  bedtools intersect -a A.bed -b B.bed
  chr1  15   20
  
.. plot::

    a = """chr1	10	20\nchr1	30	40"""
    b = """chr1	15	20"""

    title = "bedtools intersect -a A.bed -b B.bed"
    from matplotlib.pyplot import show
    from pyplots.plotter import plot_a_b_tool
    plot_a_b_tool(a, b, 'intersect', title, 'A.bed', 'B.bed')
    show()

=============================================
5.1.3 (-wa) Reporting the original A feature 
=============================================
Instead, one can force **bedtools intersect** to report the *original* **"A"** feature when an overlap is found. As
shown below, the entire "A" feature is reported, not just the portion that overlaps with the "B" feature.

For example:
::
  cat A.bed
  chr1  10  20
  chr1  30   40

  cat B.bed
  chr1  15  20

  bedtools intersect -a A.bed -b B.bed -wa
  chr1  10   20

.. plot::

    a = """chr1	10	20\nchr1	30	40"""
    b = """chr1	15	20"""

    title = "bedtools intersect -a A.bed -b B.bed -wa"
    from matplotlib.pyplot import show
    from pyplots.plotter import plot_a_b_tool
    plot_a_b_tool(a, b, 'intersect', title, 'A.bed', 'B.bed', wa=True)
    show()


=============================================
5.1.4 (-wb) Reporting the original B feature 
=============================================
Similarly, one can force **bedtools intersect** to report the *original* **"B"** feature when an overlap is found. If
just -wb is used, the overlapping portion of A will be reported followed by the *original* **"B"**. If both -wa
and -wb are used, the *originals* of both **"A"** and **"B"** will be reported.

For example (-wb alone):
::
For example:
::
  cat A.bed
  chr1  10  20
  chr1  30  40

  cat B.bed
  chr1  15   20

  bedtools intersect -a A.bed -b B.bed -wb
  chr1  15  20  chr 15  20
  

Now -wa and -wb:
::
  cat A.bed
  chr1  10  20
  chr1  30  40

  cat B.bed
  chr1  15   20

  bedtools intersect -a A.bed -b B.bed -wa -wb
  chr1  10  20  chr 15  20

=======================================================================
5.1.5 (-u) Reporting the presence of *at least one* overlapping feature 
=======================================================================
Frequently a feature in "A" will overlap with multiple features in "B". By default, **bedtools intersect** will
report each overlap as a separate output line. However, one may want to simply know that there is at
least one overlap (or none). When one uses the -u option, "A" features that overlap with one or more
"B" features are reported once. Those that overlap with no "B" features are not reported at all.


For example (*without* -u):
::
  cat A.bed
  chr1  10  20
  chr1  30  40

  cat B.bed
  chr1  15   20
  chr1  18   25
  
  bedtools intersect -a A.bed -b B.bed -wb
  chr1  10  20  chr 15  20
  chr1  10  20  chr 18  25
  
For example (*with* -u):
::
    cat A.bed
    chr1  10  20
    chr1  30  40

    cat B.bed
    chr1  15   20
    chr1  18   25

    bedtools intersect -a A.bed -b B.bed -u
    chr1  10  20

=======================================================================
5.1.6 (-c) Reporting the number of overlapping features 
=======================================================================
The -c option reports a column after each "A" feature indicating the *number* (0 or more) of overlapping
features found in "B". Therefore, *each feature in A is reported once*.

For example:
::
    cat A.bed
    chr1    10    20
    chr1    30    40

    cat B.bed
    chr1    15  20
    chr1    18  25

    bedtools intersect -a A.bed -b B.bed -u
    chr1    10    20    2
    chr1    30    40    0


=======================================================================
5.1.7 (-v) Reporting the absence of any overlapping features 
=======================================================================
There will likely be cases where you'd like to know which "A" features do not overlap with any of the
"B" features. Perhaps you'd like to know which SNPs don't overlap with any gene annotations. The -v
(an homage to "grep -v") option will only report those "A" features that have no overlaps in "B".

For example:
::
    cat A.bed
    chr1  10  20
    chr1  30  40

    cat B.bed
    chr1  15  20

    bedtools intersect -a A.bed -b B.bed -v
    chr1  30   40

.. plot::

    a = """chr1	10	20\nchr1	30	40"""
    b = """chr1	15	20"""

    title = "bedtools intersect -a A -b B -v"
    from matplotlib.pyplot import show
    from pyplots.plotter import plot_a_b_tool
    plot_a_b_tool(a, b, 'intersect', title, 'A.bed', 'B.bed', v=True)
    show()


=======================================================================
5.1.8 (-f) Requiring a minimal overlap fraction 
=======================================================================
By default, **bedtools intersect** will report an overlap between A and B so long as there is at least one base
pair is overlapping. Yet sometimes you may want to restrict reported overlaps between A and B to cases
where the feature in B overlaps at least X% (e.g. 50%) of the A feature. The -f option does exactly
this.

For example (note that the second B entry is not reported):
::
  cat A.bed
  chr1 100 200
  
  cat B.bed
  chr1 130 201
  chr1 180 220
  
  bedtools intersect -a A.bed -b B.bed -f 0.50 -wa -wb
  chr1 100 200 chr1 130 201

==========================================================================
5.1.9 (-r, combined with -f)Requiring reciprocal minimal overlap fraction 
==========================================================================
Similarly, you may want to require that a minimal fraction of both the A and the B features is
overlapped. For example, if feature A is 1kb and feature B is 1Mb, you might not want to report the
overlap as feature A can overlap at most 1% of feature B. If one set -f to say, 0.02, and one also
enable the -r (reciprocal overlap fraction required), this overlap would not be reported.

For example (note that the second B entry is not reported):
::
  cat A.bed
  chr1 100 200
  
  cat B.bed
  chr1 130 201
  chr1 130 200000
  
  bedtools intersect -a A.bed -b B.bed -f 0.50 -r -wa -wb
  chr1 100 200 chr1 130 201

==========================================================================
5.1.10 (-s)Enforcing "strandedness" 
==========================================================================
By default, **bedtools intersect** will report overlaps between features even if the features are on opposite
strands. However, if strand information is present in both BED files and the "-s" option is used, overlaps
will only be reported when features are on the same strand.

For example (note that the second B entry is not reported):
::
  cat A.bed
  chr1 100 200 a1 100 +
  
  cat B.bed
  chr1 130 201 b1 100 -
  chr1 130 201 b2 100 +
  
  bedtools intersect -a A.bed -b B.bed -wa -wb -s
  chr1 100 200 a1 100 + chr1 130 201 b2 100 +
  
  
==========================================================================
5.1.11 (-abam)Default behavior when using BAM input 
==========================================================================
When comparing alignments in BAM format (**-abam**) to features in BED format (**-b**), **bedtools intersect**
will, **by default**, write the output in BAM format. That is, each alignment in the BAM file that meets
the user's criteria will be written (to standard output) in BAM format. This serves as a mechanism to
create subsets of BAM alignments are of biological interest, etc. Note that only the mate in the BAM
alignment is compared to the BED file. Thus, if only one end of a paired-end sequence overlaps with a
feature in B, then that end will be written to the BAM output. By contrast, the other mate for the
pair will not be written. One should use **pairToBed(Section 5.2)** if one wants each BAM alignment
for a pair to be written to BAM output.

For example:
::
  bedtools intersect -abam reads.unsorted.bam -b simreps.bed | samtools view - | head -3
  
  BERTHA_0001:3:1:15:1362#0 99 chr4 9236904 0 50M = 9242033 5 1 7 9
  AGACGTTAACTTTACACACCTCTGCCAAGGTCCTCATCCTTGTATTGAAG W c T U ] b \ g c e g X g f c b f c c b d d g g V Y P W W _
  \c`dcdabdfW^a^gggfgd XT:A:R NM:i:0 SM:i:0 AM:i:0 X0:i:19 X1:i:2 XM:i:0 XO:i:0 XG:i:0 MD:Z:50
  BERTHA _0001:3:1:16:994#0 83 chr6 114221672 37 25S6M1I11M7S =
  114216196 -5493 G A A A G G C C A G A G T A T A G A A T A A A C A C A A C A A T G T C C A A G G T A C A C T G T T A
  gffeaaddddggggggedgcgeggdegggggffcgggggggegdfggfgf XT:A:M NM:i:3 SM:i:37 AM:i:37 XM:i:2 X O : i :
  1 XG:i:1 MD:Z:6A6T3
  BERTHA _0001:3:1:16:594#0 147 chr8 43835330 0 50M =
  43830893 -4487 CTTTGGGAGGGCTTTGTAGCCTATCTGGAAAAAGGAAATATCTTCCCATG U
  \e^bgeTdg_Kgcg`ggeggg_gggggggggddgdggVg\gWdfgfgff XT:A:R NM:i:2 SM:i:0 AM:i:0 X0:i:10 X1:i:7 X M : i :
  2 XO:i:0 XG:i:0 MD:Z:1A2T45
  

==========================================================================
5.1.12 (-bed)Output BED format when using BAM input 
==========================================================================
When comparing alignments in BAM format (**-abam**) to features in BED format (**-b**), **bedtools intersect**
will **optionally** write the output in BED format. That is, each alignment in the BAM file is converted
to a 6 column BED feature and if overlaps are found (or not) based on the user's criteria, the BAM
alignment will be reported in BED format. The BED "name" field is comprised of the RNAME field in
the BAM alignment. If mate information is available, the mate (e.g., "/1" or "/2") field will be
appended to the name. The "score" field is the mapping quality score from the BAM alignment.

For example:
::
  bedtools intersect -abam reads.unsorted.bam -b simreps.bed -bed | head -20
  
  chr4  9236903   9236953   BERTHA_0001:3:1:15:1362#0/1  0   +
  chr6  114221671 114221721 BERTHA_0001:3:1:16:994#0/1   37  -
  chr8  43835329  43835379  BERTHA_0001:3:1:16:594#0/2   0   -
  chr4  49110668  49110718  BERTHA_0001:3:1:31:487#0/1   23  +
  chr19 27732052  27732102  BERTHA_0001:3:1:32:890#0/2   46  +
  chr19 27732012  27732062  BERTHA_0001:3:1:45:1135#0/1  37  +
  chr10 117494252 117494302 BERTHA_0001:3:1:68:627#0/1   37  -
  chr19 27731966  27732016  BERTHA_0001:3:1:83:931#0/2   9   +
  chr8  48660075  48660125  BERTHA_0001:3:1:86:608#0/2   37  -
  chr9  34986400  34986450  BERTHA_0001:3:1:113:183#0/2  37  -
  chr10 42372771  42372821  BERTHA_0001:3:1:128:1932#0/1 3   -
  chr19 27731954  27732004  BERTHA_0001:3:1:130:1402#0/2 0   +
  chr10 42357337  42357387  BERTHA_0001:3:1:137:868#0/2  9   +
  chr1  159720631 159720681 BERTHA_0001:3:1:147:380#0/2  37  -
  chrX  58230155  58230205  BERTHA_0001:3:1:151:656#0/2  37  -
  chr5  142612746 142612796 BERTHA_0001:3:1:152:1893#0/1 37  -
  chr9  71795659  71795709  BERTHA_0001:3:1:177:387#0/1  37  +
  chr1  106240854 106240904 BERTHA_0001:3:1:194:928#0/1  37  -
  chr4  74128456  74128506  BERTHA_0001:3:1:221:724#0/1  37  -
  chr8  42606164  42606214  BERTHA_0001:3:1:244:962#0/1  37  +
  
==================================================================================
5.1.13 (-split)Reporting overlaps with spliced alignments or blocked BED features 
==================================================================================
As described in section 1.3.19, bedtools intersect will, by default, screen for overlaps against the entire span
of a spliced/split BAM alignment or blocked BED12 feature. When dealing with RNA-seq reads, for
example, one typically wants to only screen for overlaps for the portions of the reads that come from
exons (and ignore the interstitial intron sequence). The **-split** command allows for such overlaps to be
performed.

For example, the diagram below illustrates the *default* behavior. The blue dots represent the "split/
spliced" portion of the alignment (i.e., CIGAR "N" operation). In this case, the two exon annotations
are reported as overlapping with the "split" BAM alignment, but in addition, a third feature that
overlaps the "split" portion of the alignment is also reported.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Exons       ---------------                                       ----------
  
  BED/BAM  A     ************.......................................****
  
  BED File B  ^^^^^^^^^^^^^^^                     ^^^^^^^^          ^^^^^^^^^^
  
  Result      ===============                     ========          ==========

  
In contrast, when using the **-split** option, only the exon overlaps are reported.
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Exons       ---------------                                       ----------
  
  BED/BAM  A     ************.......................................****
  
  BED File B  ^^^^^^^^^^^^^^^                     ^^^^^^^^          ^^^^^^^^^^
  
  Result      ===============                                       ==========