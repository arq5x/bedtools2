==============================
          BEDTools         
==============================

Created by Aaron Quinlan Spring 2009.

Copyright 2009,2010,2011 Aaron Quinlan. All rights reserved.

Stable releases: http://code.google.com/p/bedtools

Repository:      https://github.com/arq5x/bedtools

Released under GNU public license version 2 (GPL v2).


Summary
-------
BEDTools is a collection of utilities for comparing, summarizing, and 
intersecting genomic features in BED, GTF/GFF, VCF and BAM formats. 


Manual
------
See the extensive PDF manual included at: http://code.google.com/p/bedtools/downloads/detail?name=BEDTools-User-Manual.v4.pdf.

This manual covers many common usage examples.  There are also examples available at:
http://code.google.com/p/bedtools/wiki/Usage
http://code.google.com/p/bedtools/wiki/UsageAdvanced

Installation
------------
Git
...
git clone git://github.com/arq5x/bedtools.git

Download tarball - that big gray button on the upper right.
...........................................................
1. Unpack the source downloaded tarball.
2. cd into the expanded folder.
3. Type "make" and hit enter.
4. If you encountered no errors, then bedtools should now be in bin/
  If not, try to troubleshoot then email me: aaronquinlan at gmail dot com
5. Run the test suite with: "make test"
6. Copy the files in bin/ to ~/bin or if you have the privileges, to /usr/local/bin.  Make sure that the directory to which you copy the tools is in your $PATH
7. Use the tools.


