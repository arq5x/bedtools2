###########
Quick start
###########

================
Install BEDTools
================
::

  curl http://bedtools.googlecode.com/files/BEDTools.<version>.tar.gz > BEDTools.tar.gz
  tar -zxvf BEDTools.tar.gz
  cd BEDTools
  make clean
  make all
  sudo cp bin/* /usr/local/bin/

===============
Use BEDTools
===============
Below are examples of typical BEDTools usage. **Additional usage examples are described in
section 6 of this manual.** Using the "-h" option with any BEDTools will report a list of all command
line options.

A. Report the base-pair overlap between the features in two BED files.
::

  intersectBed -a reads.bed -b genes.bed

B. Report those entries in A that overlap NO entries in B. Like "grep -v"
::

  intersectBed -a reads.bed -b genes.bed ¨Cv

C. Read BED A from stdin. Useful for stringing together commands. For example, find genes that overlap LINEs
but not SINEs.
::

  intersectBed -a genes.bed -b LINES.bed | intersectBed -a stdin -b SINEs.bed ¨Cv

D. Find the closest ALU to each gene.
::

  closestBed -a genes.bed -b ALUs.bed
  
E. Merge overlapping repetitive elements into a single entry, returning the number of entries merged.
::

  mergeBed -i repeatMasker.bed -n
  
F. Merge nearby repetitive elements into a single entry, so long as they are within 1000 bp of one another.
::

  mergeBed -i repeatMasker.bed -d 1000
  
  




    