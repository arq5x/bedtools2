###############
5.8 mergeBed
###############
**mergeBed** combines overlapping or "book-ended" (that is, one base pair away) features in a feature file
into a single feature which spans all of the combined features.

==========================================================================
5.8.1 Usage and option summary 
==========================================================================
Usage:
::
  mergeBed [OPTIONS] -i <BED/GFF/VCF>
  
===========================      ===============================================================================================================================================================================================================
Option                           Description
===========================      ===============================================================================================================================================================================================================
**-s**				             Force strandedness. That is, only merge features that are the same strand. *By default, this is disabled*.
**-n**					         Report the number of BED entries that were merged. *1 is reported if no merging occurred*.
**-d**                           Maximum distance between features allowed for features to be merged. *Default is 0. That is, overlapping and/or book-ended features are merged*.
**-nms**                         Report the names of the merged features separated by semicolons.
===========================      ===============================================================================================================================================================================================================





==========================================================================
5.8.2 Default behavior
==========================================================================
Figure:
::
  Chromosome  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  BED FILE       *************   ***************   **********************
                         ********
  
  Result         ===============================   ======================
  
  
  
For example:
::
  cat A.bed
  chr1  100  200
  chr1  180  250
  chr1  250  500
  chr1  501  1000

  mergeBed -i A.bed
  chr1  100  500
  chr1  501  1000
  
  
  
  
  

==========================================================================
5.8.3 (-s)Enforcing "strandedness" 
==========================================================================
This option behaves the same as the -s option for intersectBed while scanning for features that should
be merged. Only features on the same strand will be merged. See the discussion in the intersectBed
section for details.

==========================================================================
5.8.4 (-n)Reporting the number of features that were merged 
==========================================================================
The -n option will report the number of features that were combined from the original file in order to
make the newly merged feature. If a feature in the original file was not merged with any other features,
a "1" is reported.

For example:
::
  cat A.bed
  chr1  100  200
  chr1  180  250
  chr1  250  500
  chr1  501  1000
  
  mergeBed -i A.bed -n
  chr1  100  500  3
  chr1  501  1000 1


==========================================================================
5.8.5 (-d)Controlling how close two features must be in order to merge 
==========================================================================
By default, only overlapping or book-ended features are combined into a new feature. However, one can
force mergeBed to combine more distant features with the -d option. For example, were one to set -d to
1000, any features that overlap or are within 1000 base pairs of one another will be combined.

For example:
::
  cat A.bed
  chr1  100  200
  chr1  501  1000
  
  mergeBed -i A.bed
  chr1  100  200
  chr1  501  1000

  mergeBed -i A.bed -d 1000
  chr1  100  200  1000

==========================================================================
5.8.6 (-nms)Reporting the names of the features that were merged 
==========================================================================
Occasionally, one might like to know that names of the features that were merged into a new feature.
The -nms option will add an extra column to the mergeBed output which lists (separated by
semicolons) the names of the merged features.

For example:
::
  cat A.bed
  chr1  100  200  A1
  chr1  150  300  A2
  chr1  250  500  A3
 
  mergeBed -i A.bed -nms
  chr1  100  500  A1;A2;A3
