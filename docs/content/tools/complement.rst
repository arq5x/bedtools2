.. _complement:

###############
*complement*
###############

|

.. image:: ../images/tool-glyphs/complement-glyph.png 
    :width: 600pt 

|


``bedtools complement`` returns all intervals in a genome that are **not**
covered by at least one interval in the input BED/GFF/VCF file.

    
.. seealso::

    :doc:`../tools/merge`
    

==========================================================================
Usage and option summary
==========================================================================
**Usage**:
::

  bedtools complement -i <BED/GFF/VCF> -g <GENOME>

**(or)**:
::

  complementBed -i <BED/GFF/VCF> -g <GENOME>


==========================================================================
Default behavior
==========================================================================
By default, ``bedtools complement`` returns all genomic intervals that are not
covered by at least one record from the input file.

.. code-block:: bash

  $ cat A.bed
  chr1  100  200
  chr1  400  500
  chr1  500  800

  $ cat my.genome
  chr1  1000
  chr2  800
  
  $ bedtools complement -i A.bed -g my.genome
  chr1  0   100
  chr1  200 400
  chr1  800 1000
  chr2  0   800



