###############
5.13 shuffleBed
###############
**shuffleBed** will randomly permute the genomic locations of a fearure file among a genome defined in a
genome file. One can also provide an "exclusions" BED/GFF/VCF file that lists regions where you do
not want the permuted features to be placed. For example, one might want to prevent features from
being placed in known genome gaps. **shuffleBed** is useful as a *null* basis against which to test the
significance of associations of one feature with another.



==========================================================================
5.13.1 Usage and option summary
==========================================================================
Usage:
::
  shuffleBed [OPTIONS] -i <BED/GFF/VCF> -g <GENOME>
  
===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-excl**				         A BED file of coordinates in which features from -i should *not* be placed (e.g., genome gaps).							 
**-chrom**					     Keep features in -i on the same chromosome. Solely permute their location on the chromosome. *By default, both the chromosome and position are randomly chosen*.
**-seed**                        Supply an integer seed for the shuffling. This will allow feature shuffling experiments to be recreated exactly as the seed for the pseudo-random number generation will be constant. *By default, the seed is chosen automatically*.
===========================      ===============================================================================================================================================================================================================




==========================================================================
5.13.2 Default behavior
==========================================================================
By default, **shuffleBed** will reposition each feature in the input BED file on a random chromosome at a
random position. The size and strand of each feature are preserved.

For example:
::
  cat A.bed
  chr1  0  100  a1  1  +
  chr1  0  1000 a2  2  -

  cat my.genome
  chr1  10000
  chr2  8000
  chr3  5000
  chr4  2000

  shuffleBed -i A.bed -g my.genome
  chr4  1498  1598  a1  1  +
  chr3  2156  3156  a2  2  -





==========================================================================
5.13.3 (-chrom)Requiring that features be shuffled on the same chromosome 
==========================================================================
The "**-chrom**" option behaves the same as the default behavior except that features are randomly
placed on the same chromosome as defined in the BED file.

For example:
::
  cat A.bed
  chr1  0  100  a1  1  +
  chr1  0  1000 a2  2  -

  cat my.genome
  chr1  10000
  chr2  8000
  chr3  5000
  chr4  2000

  shuffleBed -i A.bed -g my.genome -chrom
  chr1  9560  9660  a1  1  +
  chr1  7258  8258  a2  2  -

  
  
  
==========================================================================
5.13.4 Excluding certain genome regions from shuffleBed
==========================================================================
One may want to prevent BED features from being placed in certain regions of the genome. For
example, one may want to exclude genome gaps from permutation experiment. The "**-excl**" option
defines a BED file of regions that should be excluded. **shuffleBed** will attempt to permute the
locations of all features while adhering to the exclusion rules. However it will stop looking for an
appropriate location if it cannot find a valid spot for a feature after 1,000,000 tries.

For example (*note that the exclude file excludes all but 100 base pairs of the chromosome*):
::
  cat A.bed
  chr1  0  100   a1  1  +
  chr1  0  1000  a2  2  -

  cat my.genome
  chr1  10000

  cat exclude.bed
  chr1  100  10000

  shuffleBed -i A.bed -g my.genome -excl exclude.bed
  chr1  0  100  a1  1  +
  Error, line 2: tried 1000000 potential loci for entry, but could not avoid excluded
  regions. Ignoring entry and moving on.
  

For example (*now the exclusion file only excludes the first 100 bases of the chromosome*):
::
  cat A.bed
  chr1  0  100  a1  1  +
  chr1  0  1000 a2  2  -

  cat my.genome
  chr1  10000

  cat exclude.bed
  chr1  0  100

  shuffleBed -i A.bed -g my.genome -excl exclude.bed
  chr1  147  247  a1  1  +
  chr1  2441 3441 a2  2  -


==========================================================================
5.13.5 Defining a "seed" for the random replacement.
==========================================================================
**shuffleBed** uses a pseudo-random number generator to permute the locations of BED features.
Therefore, each run should produce a different result. This can be problematic if one wants to exactly
recreate an experiment. By using the "**-seed**" option, one can supply a custom integer seed for
**shuffleBed**. In turn, each execution of **shuffleBed** with the same seed and input files should produce
identical results.

For example (*note that the exclude file below excludes all but 100 base pairs of the chromosome*):
::
  cat A.bed
  chr1 0 100 a1 1 +
  chr1 0 1000 a2 2 -

  cat my.genome
  chr1 10000

  shuffleBed -i A.bed -g my.genome -seed 927442958
  chr1 6177 6277 a1 1 +
  chr1 8119 9119 a2 2 -

  shuffleBed -i A.bed -g my.genome -seed 927442958
  chr1 6177 6277 a1 1 +
  chr1 8119 9119 a2 2 -
  
  . . .
  
  shuffleBed -i A.bed -g my.genome -seed 927442958
  chr1 6177 6277 a1 1 +
  chr1 8119 9119 a2 2 -
