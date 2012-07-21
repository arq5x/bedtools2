###############
5.14 slopBed
###############
**slopBed** will increase the size of each feature in a feature file be a user-defined number of bases. While
something like this could be done with an "**awk '{OFS="\t" print $1,$2-<slop>,$3+<slop>}'**",
**slopBed** will restrict the resizing to the size of the chromosome (i.e. no start < 0 and no end >
chromosome size).


==========================================================================
5.14.1 Usage and option summary
==========================================================================
Usage:
::
  slopBed [OPTIONS] -i <BED/GFF/VCF> -g <GENOME> [-b or (-l and -r)]
  
===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-b**				             Increase the BED/GFF/VCF entry by the same number base pairs in each direction. *Integer*.							 
**-l**					         The number of base pairs to subtract from the start coordinate. *Integer*.
**-r**                           The number of base pairs to add to the end coordinate. *Integer*.
**-s**                           Define -l and -r based on strand. For example. if used, -l 500 for a negative-stranded feature, it will add 500 bp to the *end* coordinate.
===========================      ===============================================================================================================================================================================================================



==========================================================================
5.14.2 Default behavior
==========================================================================
By default, **slopBed** will either add a fixed number of bases in each direction (**-b**) or an asymmetric
number of bases in each direction (**-l** and **-r**).

For example:
::
  cat A.bed
  chr1 5 100
  chr1 800 980

  cat my.genome
  chr1 1000

  slopBed -i A.bed -g my.genome -b 5
  chr1 0 105
  chr1 795 985

  slopBed -i A.bed -g my.genome -l 2 -r 3
  chr1 3 103
  chr1 798 983
  

However, if the requested number of bases exceeds the boundaries of the chromosome, **slopBed** will
"clip" the feature accordingly.
::
  cat A.bed
  chr1  5   100
  chr1  800 980

  cat my.genome
  chr1  1000

  slopBed -i A.bed -g my.genome -b 5000
  chr1  0   1000
  chr1  0   1000

  
  
==========================================================================
5.14.3 Resizing features according to strand
==========================================================================
**slopBed** will optionally increase the size of a feature based on strand.

For example:
::
  cat A.bed
  chr1 100 200 a1 1 +
  chr1 100 200 a2 2 -

  cat my.genome
  chr1 1000

  slopBed  -i A.bed -g my.genome -l 50 -r 80 -s
  chr1 50  280 a1 1 +
  chr1 20  250 a2 2 -
