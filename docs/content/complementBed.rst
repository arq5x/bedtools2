###############
5.17 complementBed
###############
**complementBed** returns the intervals in a genome that are not by the features in a feature file. An
example usage of this tool would be to return the intervals of the genome that are not annotated as a
repeat.


==========================================================================
5.17.1 Usage and option summary
==========================================================================
Usage:
::
  complementBed [OPTIONS] -i <BED/GFF/VCF> -g <GENOME>

**No additional options.**




==========================================================================
5.17.2 Default behavior
==========================================================================
Figure:
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BED FILE A     *************   ***************     ******************              
  
  Result      ===             ===               =====                  =======


For example:
::
  cat A.bed
  chr1  100  200
  chr1  400  500
  chr1  500  800

  cat my.genome
  chr1  1000

  complementBed -i A.bed -g my.genome
  chr1  0    100
  chr1  200  400
  chr1  800  1000


