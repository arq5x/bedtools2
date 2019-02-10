############
FAQ
############


====================
Installation issues
====================

--------------------------------------------------
Why am I getting all of these zlib errors?
--------------------------------------------------

On certain operating systems (especially free Linux distributions) the complete
zlib libraries are not installed.  Bedtools depends upon zlib in order to 
decompress gzipped files.  

.. code-block:: bash

    - Building main bedtools binary.
    obj/gzstream.o: In function gzstreambuf::open(char const*, int):
    gzstream.C:(.text+0x2a5): undefined reference to gzopen64'
    collect2: ld returned 1 exit status
    make: *** [all] Error 1
    
If you see an error such as the above, it suggests you need to install the
``zlib`` and ``zlib1g-dev`` libraries.  This is typically straightforward using
package managers.  For example, on Debian/Ubuntu this would be:

.. code-block:: bash
    
    apt-get install zlib
    apt-get install zlib1g-dev

and on Fedora/Centos this would be:

.. code-block:: bash
    
    yum install zlib
    yum install zlib1g-dev

--------------------------------------------------
Compiling with a specific zlib library
--------------------------------------------------
If you need to override the location of the zlib library used to compile bedtools, you can run `make` and specify the `LIBS` argument. For example:

.. code-block:: bash

    make LIBS='/PATH/TO/ZLIB/lib/libz.a'

====================
General questions
====================


--------------------------------------------------
How do I know what version of bedtools I am using?
--------------------------------------------------

Use the --version option.

.. code-block:: bash

    $ bedtools --version
    bedtools v2.17.0


--------------------------------------------------
How do I bring up the help/usage menu?
--------------------------------------------------

To receive a high level list of available tools in bedtools, use ```-h``:

.. code-block:: bash

    $ bedtools -h
    bedtools: flexible tools for genome arithmetic and DNA sequence analysis.
    usage:    bedtools <subcommand> [options]
    
    The bedtools sub-commands include:
    
    [ Genome arithmetic ]
        intersect     Find overlapping intervals in various ways.
        window        Find overlapping intervals within a window around an interval.
        closest       Find the closest, potentially non-overlapping interval.
        coverage      Compute the coverage over defined intervals.
        map           Apply a function to a column for each overlapping interval.
        genomecov     Compute the coverage over an entire genome.
        merge         Combine overlapping/nearby intervals into a single interval.
        cluster       Cluster (but don't merge) overlapping/nearby intervals.
        complement    Extract intervals _not_ represented by an interval file.
    ...

To display the help for a specific tool (e.g., ``bedtools shuffle``), use:

.. code-block:: bash

    $ bedtools merge -h
    
    Tool:    bedtools merge (aka mergeBed)
    Version: v2.17.0
    Summary: Merges overlapping BED/GFF/VCF entries into a single interval.
    
    Usage:   bedtools merge [OPTIONS] -i <bed/gff/vcf>
    
    Options: 
    	-s	Force strandedness.  That is, only merge features
    		that are the same strand.
    		- By default, merging is done without respect to strand.
    
    	-n	Report the number of BED entries that were merged.
    		- Note: "1" is reported if no merging occurred.



            

====================
Issues with output
====================

------------------------------------------------------------------------
I *know* there are overlaps, but none are reported. What might be wrong?
------------------------------------------------------------------------

There are two common causes of this problem.  The first cause is non-obvious 
differences in the way chromosomes are named in files being compared.  
For example, "1" is not the same as "chr1" just as "   chr1" is not the same 
as "chr1".  Secondly, users often copy files from a Windows machine to a UNIX 
machine.  This causes issues because Windows uses two bytes to represent
the end of a line (``\r\n``) whereas the UNIX convention uses a single byte
(``\n``).  If your files don't conform to the UNIX convention, you will have 
problems.  One can convert files from Windows to UNIX with
the following command:

.. code-block:: bash

   perl -i -p -e 's/\r\n/\n/g;' file.windows > file.unix



====================
Installation issues
====================


---------------------------------------------------------------------------
Bedtools compilation fails with errors related to zlib.  How do I fix this?
---------------------------------------------------------------------------

Some systems, especially Ubuntu, do not come pre-installed with up to date
versions of the zlib compression utilities that tools such as `bedtools` and
`samtools` depend upon. This can cause compilation errors when you try to 
compile `bedtools`.  Errors include:

.. code-block:: bash

    ../utils//gzstream/gzstream.h:50: error: ‘gzFile’ does not name a type 
    

or

.. code-block:: bash

    fatal error: zlib.h: No such file or directory  

This indicates that you need to install the zlib libraries on your system, which
turns out to not be too difficult through the use of package installers.  For
example, on Ubuntu, you'd want to run:

.. code-block:: bash

    sudo apt-get install zlib1g-dev
    sudo apt-get install zlib


