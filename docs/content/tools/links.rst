.. _links:

###############
*links*
###############
Creates an HTML file with links to an instance of the UCSC Genome Browser for all features /
intervals in a file. This is useful for cases when one wants to manually inspect through a large set of
annotations or features.

==========================================================================
Usage and option summary
==========================================================================
Usage:

::

  linksBed [OPTIONS] -i <BED/GFF/VCF> > <HTML file>
  
===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-base**				         The "basename" for the UCSC browser. *Default: http://genome.ucsc.edu*					 
**-org**					     The organism (e.g. mouse, human). *Default: human*
**-db**                          The genome build. *Default: hg18*
===========================      ===============================================================================================================================================================================================================




==========================================================================
Default behavior
==========================================================================
By default, **linksBed** creates links to the public UCSC Genome Browser.

For example:

::

  head genes.bed
  chr21 9928613  10012791  uc002yip.1 0  -
  chr21 9928613  10012791  uc002yiq.1 0  -
  chr21 9928613  10012791  uc002yir.1 0  -
  chr21 9928613  10012791  uc010gkv.1 0  -
  chr21 9928613  10061300  uc002yis.1 0  -
  chr21 10042683 10120796  uc002yit.1 0  -
  chr21 10042683 10120808  uc002yiu.1 0  -
  chr21 10079666 10120808  uc002yiv.1 0  -
  chr21 10080031 10081687  uc002yiw.1 0  -
  chr21 10081660 10120796  uc002yix.2 0  -

  linksBed -i genes.bed > genes.html
  
When genes.html is opened in a web browser, one should see something like the following, where each
link on the page is built from the features in genes.bed:





==========================================================================
Creating HTML links to a local UCSC Browser installation
==========================================================================
Optionally, **linksBed** will create links to a local copy of the UCSC Genome Browser.

For example:

::

  head -3 genes.bed
  chr21 9928613 10012791 uc002yip.1 0 -
  chr21 9928613 10012791 uc002yiq.1 0 -

  linksBed -i genes.bed -base http://mirror.uni.edu > genes.html
  
One can point the links to the appropriate organism and genome build as well:

::

  head -3 genes.bed
  chr21 9928613 10012791 uc002yip.1 0 -
  chr21 9928613 10012791 uc002yiq.1 0 -

  linksBed -i genes.bed -base http://mirror.uni.edu -org mouse -db mm9 > genes.html

