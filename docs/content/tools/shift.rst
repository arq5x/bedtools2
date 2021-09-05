.. _shift:

###############
*shift*
###############

|

.. image:: ../images/tool-glyphs/shift-glyph.png 
    :width: 600pt 

|

``bedtools shift`` will move each feature in a feature file by a 
user-defined number of bases. While something like this could be done with an 
``awk '{OFS="\t" print $1,$2+<shift>,$3+<shift>}'``,
``bedtools shift`` will restrict the resizing to the size of the chromosome 
(i.e. no features before 0 or past the chromosome end).

.. note::

    In order to prevent the extension of intervals beyond chromosome boundaries,
    ``bedtools shift`` requires a *genome* file defining the length of each 
    chromosome or contig.

.. seealso::

    :doc:`../tools/slop`

==========================================================================
Usage and option summary
==========================================================================
**Usage**:
::

  bedtools shift [OPTIONS] -i <BED/GFF/VCF> -g <GENOME> [-s or (-m and -p)]

**(or):**
::

  shiftBed [OPTIONS] -i <BED/GFF/VCF> -g <GENOME> [-s or (-m and -p)]
    
===========================      ===============================================================================================================================================================================================================
 Option                           Description
===========================      ===============================================================================================================================================================================================================
**-s**                           Shift the BED/GFF/VCF entry -s base pairs. *Integer* or *Float* (e.g. 0.1) if used with -pct
**-m**                           Shift entries on the - strand -m base pairs. *Integer* or *Float* (e.g. 0.1) if used with -pct
**-p**                           Shift entries on the + strand -p base pairs. *Integer* or *Float* (e.g. 0.1) if used with -pct
**-pct**                         Define -s, -m and -p as a fraction of the feature's length. E.g. if used on a 1000bp feature, -s 0.50, will shift the feature 500 bp \"upstream\".  Default = false.
**-header**                      Print the header from the input file prior to results.
===========================      ===============================================================================================================================================================================================================



==========================================================================
Default behavior
==========================================================================
By default, ``bedtools shift`` will shift features a fixed number of bases. This shift can either be applied to all features (``-s``), or seperately to features on the plus (``-p``) and minus (``-m``) strands. Shifts can either be positive or negative.

.. code-block:: bash

  $ cat A.bed
  chr1 5 100 +
  chr1 800 980 -

  $ cat my.genome
  chr1 1000

  $ bedtools shift -i A.bed -g my.genome -s 5
  chr1 10 105 +
  chr1 805 985 -

  $ bedtools shift -i A.bed -g my.genome -p 2 -m -3
  chr1 7 102 +
  chr1 797 977 -
  

However, if the requested number of bases exceeds the boundaries of the 
chromosome, ``bedtools shift`` will "clip" the feature accordingly.

.. code-block:: bash

  $ cat A.bed
  chr1  5   100 +
  chr1  800 980 +

  $ cat my.genome
  chr1  1000

  $ bedtools shift -i A.bed -g my.genome -s 5000
  chr1  999   1000 +
  chr1  999   1000 +


==========================================================================
``-pct`` Shifting features by a given fraction
==========================================================================
``bedtools shift`` will shift the feature by a user-specific fraction of the
feature length. Hence, 0.5 will add half the length of the feature to the start and end coordinates. 

For example:

.. code-block:: bash

  $ cat A.bed
  chr1 100 200 a1 1 +

  $ bedtools shift -i A.bed -g my.genome -s 0.5 -pct
  chr1 150  250 a1 1 +


==========================================================================
``-header`` Print the header for the A file before reporting results.
==========================================================================
By default, if your A file has a header, it is ignored when reporting results.
This option will instead tell bedtools to first print the header for the
A file prior to reporting results.

