###############
General usage
###############

=======================
Supported file formats
=======================

----------------------
BED format
----------------------
As described on the UCSC Genome Browser website (see link below), the browser extensible data (BED) format is a concise and
flexible way to represent genomic features and annotations. The BED format description supports up to
12 columns, but only the first 3 are required for the UCSC browser, the Galaxy browser and for
bedtools. bedtools allows one to use the "BED12" format (that is, all 12 fields listed below).
However, only intersectBed, coverageBed, genomeCoverageBed, and bamToBed will obey the BED12
"blocks" when computing overlaps, etc., via the **"-split"** option. For all other tools, the last six columns
are not used for any comparisons by the bedtools. Instead, they will use the entire span (start to end)
of the BED12 entry to perform any relevant feature comparisons. The last six columns will be reported
in the output of all comparisons.

The file description below is modified from: http://genome.ucsc.edu/FAQ/FAQformat#format1.

1. **chrom** - The name of the chromosome on which the genome feature exists.

  - *Any string can be used*. For example, "chr1", "III", "myChrom", "contig1112.23".
  - *This column is required*.

2. **start** - The zero-based starting position of the feature in the chromosome.

 - *The first base in a chromosome is numbered 0*.
 - *The start position in each BED feature is therefore interpreted to be 1 greater than the start position listed in the feature. For example, start=9, end=20 is interpreted to span bases 10 through 20,inclusive*.
 - *This column is required*.

3. **end** - The one-based ending position of the feature in the chromosome.

 - *The end position in each BED feature is one-based. See example above*.
 - *This column is required*.

4. **name** - Defines the name of the BED feature.

 - *Any string can be used*. For example, "LINE", "Exon3", "HWIEAS_0001:3:1:0:266#0/1", or "my_Feature".
 - *This column is optional*.

5. **score** - The UCSC definition requires that a BED score range from 0 to 1000, inclusive. However, bedtools allows any string to be stored in this field in order to allow greater flexibility in annotation features. For example, strings allow scientific notation for p-values, mean enrichment values, etc. It should be noted that this flexibility could prevent such annotations from being correctly displayed on the UCSC browser.

 - *Any string can be used*. For example, 7.31E-05 (p-value), 0.33456 (mean enrichment value), "up", "down", etc.
 - *This column is optional*.

6. **strand** - Defines the strand - either '+' or '-'.

 - *This column is optional*.

7. **thickStart** - The starting position at which the feature is drawn thickly.

 - *Allowed yet ignored by bedtools*.

8. **thickEnd** - The ending position at which the feature is drawn thickly.

 - *Allowed yet ignored by bedtools*.

9. **itemRgb** - An RGB value of the form R,G,B (e.g. 255,0,0).
 
 - *Allowed yet ignored by bedtools*.

10. **blockCount** - The number of blocks (exons) in the BED line.
 
 - *Allowed yet ignored by bedtools*.

11. **blockSizes** - A comma-separated list of the block sizes.


12. **blockStarts** - A comma-separated list of block starts.

 - *Allowed yet ignored by bedtools*.
 
 
bedtools requires that all BED input files (and input received from stdin) are **tab-delimited**. The following types of BED files are supported by bedtools:


1.  **BED3**: A BED file where each feature is described by **chrom**, **start**, and **end**.

  For example: ``chr1          11873   14409``

2.  **BED4**: A BED file where each feature is described by **chrom**, **start**, **end**, and **name**.

  For example: ``chr1  11873  14409  uc001aaa.3``

3.  **BED5**: A BED file where each feature is described by **chrom**, **start**, **end**, **name**, and **score**.
  
  For example: ``chr1 11873 14409 uc001aaa.3 0``

4.  **BED6**: A BED file where each feature is described by **chrom**, **start**, **end**, **name**, **score**, and **strand**.

  For example: ``chr1 11873 14409 uc001aaa.3 0 +``

5.  **BED12**: A BED file where each feature is described by all twelve columns listed above.

    For example: ``chr1 11873 14409 uc001aaa.3 0 + 11873 11873 0 3 354,109,1189, 0,739,1347,``

----------------------
BEDPE format
----------------------
We have defined a new file format, the browser extensible data paired-end (BEDPE) format, in order to concisely describe disjoint genome features,
such as structural variations or paired-end sequence alignments. We chose to define a new format
because the existing "blocked" BED format (a.k.a. BED12) does not allow inter-chromosomal feature
definitions. In addition, BED12 only has one strand field, which is insufficient for paired-end sequence
alignments, especially when studying structural variation.

The BEDPE format is described below. The description is modified from: http://genome.ucsc.edu/FAQ/FAQformat#format1.

1. **chrom1** - The name of the chromosome on which the **first** end of the feature exists.

 - *Any string can be used*. For example, "chr1", "III", "myChrom", "contig1112.23".
 - *This column is required*.
 - *Use "." for unknown*.

2. **start1** - The zero-based starting position of the **first** end of the feature on **chrom1**.
 
 - *The first base in a chromosome is numbered 0*.
 - *As with BED format, the start position in each BEDPE feature is therefore interpreted to be 1 greater than the start position listed in the feature. This column is required*.
 - *Use -1 for unknown*.

3. **end1** - The one-based ending position of the first end of the feature on **chrom1**.

 - *The end position in each BEDPE feature is one-based*.
 - *This column is required*.
 - *Use -1 for unknown*.

4. **chrom2** - The name of the chromosome on which the **second** end of the feature exists.

 - *Any string can be used*. For example, "chr1", "III", "myChrom", "contig1112.23".
 - *This column is required*.
 - *Use "." for unknown*.

5. **start2** - The zero-based starting position of the **second** end of the feature on **chrom2**.

 - *The first base in a chromosome is numbered 0*.
 - *As with BED format, the start position in each BEDPE feature is therefore interpreted to be 1 greater than the start position listed in the feature. This column is required*.
 - *Use -1 for unknown*.

6. **end2** - The one-based ending position of the **second** end of the feature on **chrom2**.

 - *The end position in each BEDPE feature is one-based*.
 - *This column is required*.
 - *Use -1 for unknown*.

7. **name** - Defines the name of the BEDPE feature.

 - *Any string can be used*. For example, "LINE", "Exon3", "HWIEAS_0001:3:1:0:266#0/1", or "my_Feature".
 - *This column is optional*.

8. **score** - The UCSC definition requires that a BED score range from 0 to 1000, inclusive. *However, bedtools allows any string to be stored in this field in order to allow greater flexibility in annotation features*. For example, strings allow scientific notation for p-values, mean enrichment values, etc. It should be noted that this flexibility could prevent such annotations from being correctly displayed on the UCSC browser.

 - *Any string can be used*. For example, 7.31E-05 (p-value), 0.33456 (mean enrichment value), "up", "down", etc.
 - *This column is optional*.

9. **strand1** - Defines the strand for the first end of the feature. Either '+' or '-'.

 - *This column is optional*.
 - *Use "." for unknown*.

10. **strand2** - Defines the strand for the second end of the feature. Either '+' or '-'.

 - *This column is optional*.
 - *Use "." for unknown*.

11. **Any number of additional, user-defined fields** - bedtools allows one to add as many additional fields to the normal, 10-column BEDPE format as necessary. These columns are merely "passed through" **pairToBed** and **pairToPair** and are not part of any analysis. One would use these additional columns to add extra information (e.g., edit distance for each end of an alignment, or "deletion", "inversion", etc.) to each BEDPE feature.

 - *These additional columns are optional*.

 
Entries from an typical BEDPE file:
::

  chr1  100   200   chr5  5000  5100  bedpe_example1  30   +  -
  chr9  1000  5000  chr9  3000  3800  bedpe_example2  100  +  -


Entries from a BEDPE file with two custom fields added to each record:
::

  chr1  10    20    chr5  50    60    a1     30       +    -  0  1
  chr9  30    40    chr9  80    90    a2     100      +    -  2  1



----------------------
GFF format
----------------------
The GFF format is described on the Sanger Institute's website (http://www.sanger.ac.uk/resources/software/gff/spec.html). The GFF description below is modified from the definition at this URL. All nine columns in the GFF format description are required by bedtools.

1. **seqname** - The name of the sequence (e.g. chromosome) on which the feature exists.

 - *Any string can be used*. For example, "chr1", "III", "myChrom", "contig1112.23".
 - *This column is required*.

2. **source** - The source of this feature. This field will normally be used to indicate the program making the prediction, or if it comes from public database annotation, or is experimentally verified, etc.

 - *This column is required*.

3. **feature** - The feature type name. Equivalent to BED's **name** field.

 - *Any string can be used*. For example, "exon", etc.
 - *This column is required*.

4. **start** - The one-based starting position of feature on **seqname**.
 
 - *This column is required*. 
 - *bedtools accounts for the fact the GFF uses a one-based position and BED uses a zero-based start position*.

5. **end** - The one-based ending position of feature on **seqname**.

 - *This column is required*.

6. **score** - A score assigned to the GFF feature. Like BED format, bedtools allows any string to be stored in this field in order to allow greater flexibility in annotation features. We note that this differs from the GFF definition in the interest of flexibility.

 - *This column is required*.

7. **strand** - Defines the strand. Use '+', '-' or '.'

 - *This column is required*.

8. **frame** -  The frame of the coding sequence. Use '0', '1', '2', or '.'.

 - *This column is required*.

9. **attribute** - Taken from http://www.sanger.ac.uk/resources/software/gff/spec.html: From version 2 onwards, the attribute field must have an tag value structure following the syntax used within objects in a .ace file, flattened onto one line by semicolon separators. Free text values must be quoted with double quotes. *Note: all non-printing characters in such free text value strings (e.g. newlines, tabs, control characters, etc) must be explicitly represented by their C (UNIX) style backslash-escaped representation (e.g. newlines as '\n', tabs as '\t')*. As in ACEDB, multiple values can follow a specific tag. The aim is to establish consistent use of particular tags, corresponding to an underlying implied ACEDB model if you want to think that way (but acedb is not required).

 - *This column is required*.

An entry from an example GFF file :

::

  seq1 BLASTX similarity 101 235 87.1 + 0 Target "HBA_HUMAN" 11 55 ;
  E_value 0.0003 dJ102G20 GD_mRNA coding_exon 7105 7201 . - 2 Sequence
  "dJ102G20.C1.1"
  
  
  
------------------------
*Genome* file format
------------------------
Some of the bedtools (e.g., genomeCoverageBed, complementBed, slopBed) need to know the size of
the chromosomes for the organism for which your BED files are based. When using the UCSC Genome
Browser, Ensemble, or Galaxy, you typically indicate which which species/genome build you are
working. The way you do this for bedtools is to create a "genome" file, which simply lists the names of
the chromosomes (or scaffolds, etc.) and their size (in basepairs).


Genome files must be **tab-delimited** and are structured as follows (this is an example for *C. elegans*):

::

  chrI  15072421
  chrII 15279323 
  ...
  chrX  17718854
  chrM  13794

bedtools includes pre-defined genome files for human and mouse in the **/genomes** directory included
in the bedtools distribution.

One can also create a suitable genome file by running `samtools faidx` on the appropriate
FASTA reference genome. Then use the resulting .fai file as a genome file, as bedtools will only
care about the first two columns, which define the chromosome name and length.
For example:

::

  # download GRCh38
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
  # create an index of it
  samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
  # use the .fai index as a genome file with bedtools
  bedtools complement my.grch38.bed -g GRCh38_full_analysis_set_plus_decoy_hla.fa.fai


----------------------
SAM/BAM format
----------------------
The SAM / BAM format is a powerful and widely-used format for storing sequence alignment data (see
http://samtools.sourceforge.net/ for more details). It has quickly become the standard format to which
most DNA sequence alignment programs write their output. Currently, the following bedtools
support input in BAM format: ``intersect``, ``window``, ``coverage``, ``genomecov``,
``pairtobed``, ``bamtobed``. Support for the BAM format in bedtools allows one to (to name a few):
compare sequence alignments to annotations, refine alignment datasets, screen for potential mutations
and compute aligned sequence coverage.



----------------------
VCF format
----------------------
The Variant Call Format (VCF) was conceived as part of the 1000 Genomes Project as a standardized
means to report genetic variation calls from SNP, INDEL and structural variant detection programs
(see http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf4.0 for details).
bedtools now supports the latest version of this format (i.e, Version 4.0). As a result, bedtools can
be used to compare genetic variation calls with other genomic features.
