.. _unionbedg:

###############
*unionbedg*
###############
**bedtools unionbedg** combines multiple BEDGRAPH files into a single file such that one can directly
compare coverage (and other text-values such as genotypes) across multiple sample


==========================================================================
Usage and option summary
==========================================================================
Usage:

::

  bedtools unionbedg [OPTIONS] -i FILE1 FILE2 FILE3 ... FILEn
  
===========================      ===============================================================================================================================================================================================================
 Option                           Description
 
===========================      ===============================================================================================================================================================================================================
**-header**				         Print a header line, consisting of chrom, start, end followed by the names of each input BEDGRAPH file.	 
**-names**					     A list of names (one per file) to describe each file in -i. These names will be printed in the header line.
**-empty**                       Report empty regions (i.e., start/end intervals w/o values in all files). *Requires the '-g FILE' parameter (see below)*.
**-g**                           The genome file to be used to calculate empty regions.
**-filler TEXT**                 Use TEXT when representing intervals having no value. Default is '0', but you can use 'N/A' or any other text.
**-examples**                    Show detailed usage examples.
===========================      ===============================================================================================================================================================================================================




==========================================================================
Default behavior
==========================================================================

::

  cat 1.bg
  chr1 1000 1500 10
  chr1 2000 2100 20

  cat 2.bg
  chr1 900 1600 60
  chr1 1700 2050 50

  cat 3.bg
  chr1 1980 2070 80
  chr1 2090 2100 20

  cat sizes.txt
  chr1 5000

  bedtools unionbedg -i 1.bg 2.bg 3.bg
  chr1 900  1000 0  60 0
  chr1 1000 1500 10 60 0
  chr1 1500 1600 0  60 0
  chr1 1700 1980 0  50 0
  chr1 1980 2000 0  50 80
  chr1 2000 2050 20 50 80
  chr1 2050 2070 20 0  80
  chr1 2070 2090 20 0  0
  chr1 2090 2100 20 0  20

==========================================================================
``-header`` Add a header line to the output
==========================================================================

::

  bedtools unionbedg -i 1.bg 2.bg 3.bg -header
  chrom  start  end  1  2  3
  chr1   900    1000 0  60 0
  chr1   1000   1500 10 60 0
  chr1   1500   1600 0  60 0
  chr1   1700   1980 0  50 0
  chr1   1980   2000 0  50 80
  chr1   2000   2050 20 50 80
  chr1   2050   2070 20 0  80
  chr1   2070   2090 20 0  0
  chr1   2090   2100 20 0  20


==========================================================================
``-names`` Add a header line with custom file names to the output
==========================================================================

::

  bedtools unionbedg -i 1.bg 2.bg 3.bg -header -names WT-1 WT-2 KO-1
  chrom  start  end   WT-1  WT-2  KO-1
  chr1   900    1000  0     60    0
  chr1   1000   1500  10    60    0
  chr1   1500   1600  0     60    0
  chr1   1700   1980  0     50    0
  chr1   1980   2000  0     50    80
  chr1   2000   2050  20    50    80
  chr1   2050   2070  20    0     80
  chr1   2070   2090  20    0     0
  chr1   2090   2100  20    0     20


  
  
==========================================================================
``-empty`` Include regions that have zero coverage in all BEDGRAPH files.
==========================================================================

::

  bedtools unionbedg -i 1.bg 2.bg 3.bg -empty -g sizes.txt -header
  chrom  start  end  WT-1  WT-2  KO-1
  chrom  start  end  1     2     3
  chr1   0      900  0     0     0
  chr1   900    1000 0     60    0
  chr1   1000   1500 10    60    0
  chr1   1500   1600 0     60    0
  chr1   1600   1700 0     0     0
  chr1   1700   1980 0     50    0
  chr1   1980   2000 0     50    80
  chr1   2000   2050 20    50    80
  chr1   2050   2070 20    0     80
  chr1   2070   2090 20    0     0
  chr1   2090   2100 20    0     20
  chr1   2100   5000 0     0     0


==========================================================================
``-filler`` Use a custom value for missing values.
==========================================================================

::

  bedtools unionbedg -i 1.bg 2.bg 3.bg -empty -g sizes.txt -header -filler N/A
  chrom start end  WT-1  WT-2  KO-1
  chrom start end  1     2     3
  chr1  0     900  N/A   N/A   N/A
  chr1  900   1000 N/A   60    N/A
  chr1  1000  1500 10    60    N/A
  chr1  1500  1600 N/A   60    N/A
  chr1  1600  1700 N/A   N/A   N/A
  chr1  1700  1980 N/A   50    N/A
  chr1  1980  2000 N/A   50    80
  chr1  2000  2050 20    50    80
  chr1  2050  2070 20    N/A   80
  chr1  2070  2090 20    N/A   N/A
  chr1  2090  2100 20    N/A   20
  chr1  2100  5000 N/A   N/A   N/A

  
==========================================================================
Use BEDGRAPH files with non-numeric values.
==========================================================================

::

  cat 1.snp.bg
  chr1 0 1 A/G
  chr1 5 6 C/T

  cat 2.snp.bg
  chr1 0 1 C/C
  chr1 7 8 T/T

  cat 3.snp.bg
  chr1 0 1 A/G
  chr1 5 6 C/T

  bedtools unionbedg -i 1.snp.bg 2.snp.bg 3.snp.bg -filler -/-
  chr1 0 1 A/G C/C A/G
  chr1 5 6 C/T -/- C/T
  chr1 7 8 -/- T/T -/-
