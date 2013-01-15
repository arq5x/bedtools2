############
FAQ
############


====================
Installation issues
====================

Why am I getting all of these zlib errors?
..........................................

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
