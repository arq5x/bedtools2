###############
5.7 subtractBed
###############
**subtractBed** searches for features in B that overlap A. If an overlapping feature is found in B, the
overlapping portion is removed from A and the remaining portion of A is reported. If a feature in B
overlaps all of a feature in A, the A feature will not be reported.


==========================================================================
5.7.1 Usage and option summary
==========================================================================
Usage:
::
  subtractBed [OPTIONS] -a <BED/GFF/VCF> -b <BED/GFF/VCF>
  
===========================      ===============================================================================================================================================================================================================
Option                           Description
===========================      ===============================================================================================================================================================================================================
**-f**				             Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
**-s**					         Force strandedness. That is, find the closest feature in B overlaps A on the same strand.  *By default, this is disabled*.
===========================      ===============================================================================================================================================================================================================



==========================================================================
5.7.2 Default behavior
========================================================================== 
Figure:
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BED FILE A             *************            ******
  
  BED File B         ^^^^^^^^                   ^^^^^^^^^^^
  
  Result                     =========
  
For example:
::
  cat A.bed
  chr1  100  200
  chr1  10   20

  cat B.bed
  chr1  0    30
  chr1  180  300

  subtractBed -a A.bed -b B.bed
  chr1  100  180
  
  
  
  
  

==========================================================================
5.7.3  (-f)Requiring a minimal overlap fraction before subtracting
==========================================================================
This option behaves the same as the -f option for intersectBed. In this case, subtractBed will only
subtract an overlap with B if it covers at least the fraction of A defined by -f. If an overlap is found,
but it does not meet the overlap fraction, the original A feature is reported without subtraction.

For example:
::
  cat A.bed
  chr1  100  200

  cat B.bed
  chr1  180  300

  subtractBed -a A.bed -b B.bed -f 0.10
  chr1  100  180

  subtractBed -a A.bed -b B.bed -f 0.80
  chr1  100  200




==========================================================================
5.7.4 (-s)Enforcing "strandedness" 
==========================================================================
This option behaves the same as the -s option for intersectBed while scanning for features in B that
should be subtracted from A. See the discussion in the intersectBed section for details.




