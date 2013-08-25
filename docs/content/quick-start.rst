###########
Quick start
###########

================
Install bedtools
================

.. code-block:: bash

  curl http://bedtools.googlecode.com/files/BEDTools.<version>.tar.gz > BEDTools.tar.gz
  tar -zxvf BEDTools.tar.gz
  cd BEDTools
  make
  sudo cp bin/* /usr/local/bin/


===============
Use bedtools
===============
Below are examples of typical bedtools usage. Using the "-h" option with any 
bedtools will report a list of all command line options.


Report the base-pair overlap between the features in two BED files.

.. code-block:: bash

  bedtools intersect -a reads.bed -b genes.bed


Report those entries in A that overlap NO entries in B. Like "grep -v"

.. code-block:: bash

  bedtools intersect  -a reads.bed -b genes.bed -v


Read BED A from STDIN. Useful for stringing together commands. For example, 
find genes that overlap LINEs but not SINEs.

.. code-block:: bash

  bedtools intersect -a genes.bed -b LINES.bed | \
    bedtools intersect -a stdin -b SINEs.bed -v


Find the closest ALU to each gene.

.. code-block:: bash

   bedtools closest -a genes.bed -b ALUs.bed
  

Merge overlapping repetitive elements into a single entry, returning the number of entries merged.

.. code-block:: bash

  bedtools merge -i repeatMasker.bed -n
  

Merge nearby repetitive elements into a single entry, so long as they are within 1000 bp of one another.

.. code-block:: bash

  bedtools merge -i repeatMasker.bed -d 1000
  
  




    
