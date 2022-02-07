############
Installation
############


``bedtools`` is intended to run in a "command line" environment on UNIX, LINUX
and Apple OS X operating systems. Installing ``bedtools`` involves either
downloading the source code and compiling it manually, or installing stable
release from package managers such as
`homebrew (for OS X) <http://mxcl.github.com/homebrew/>`_.




--------------------------
Installing stable releases
--------------------------

.....................................
Downloading a pre-compiled binary
.....................................
.. note::    
    1. The following approach will only work for Linux (non-OSX) systems.

Starting with release 2.28.0, we provide statically-linked binaries thast should work
right away on Linux systems. Go to the `releases <https://github.com/arq5x/bedtools2/releases>`_
page and look for the static binary named `bedtools.static.binary`. Right click on it, get the URL, then download it
with `wget` or `curl` and you should be good to go.

As an example. First, download `bedtools.static.binary` from the latest release `here <https://github.com/arq5x/bedtools2/releases>`_

.. code-block:: bash

    mv bedtools.static.binary bedtools
    chmod a+x bedtools


.....................................
Compiling from source via Github
.....................................

Stable, versioned releases of bedtools are made available on Github at the
`bedtools2 repository <https://github.com/arq5x/bedtools2/>`_ under 
the `releases <https://github.com/arq5x/bedtools2/releases>`_ tab.  
The following commands will install ``bedtools`` in a local directory on an UNIX or OS X machine.

.. note::    
    1. The bedtools Makefiles utilize the GCC compiler. One should edit the
    Makefiles accordingly if one wants to use a different compiler.

.. code-block:: bash

  $ wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz 
  $ tar -zxvf bedtools-2.29.1.tar.gz
  $ cd bedtools2
  $ make

At this point, one should copy the binaries in ./bin/ to either
``usr/local/bin/`` or some other repository for commonly used UNIX tools in
your environment. You will typically require administrator (e.g. "root" or
"sudo") privileges to copy to ``usr/local/bin/``. If in doubt, contact you
system administrator for help.

.....................................
Installing with package managers
.....................................

In addition, stable releases of ``bedtools`` are also available through package
managers such as `homebrew (for OS X) <http://mxcl.github.com/homebrew/>`_,
``apt-get`` and ``yum``.

**Fedora/Centos**. Adam Huffman has created a Red Hat package for bedtools so
that one can easily install the latest release using "yum", the Fedora
package manager. It should work with Fedora 13, 14 and EPEL5/6 (
for Centos, Scientific Linux, etc.).

.. code-block:: bash

    yum install BEDTools

**Debian/Ubuntu.** Charles Plessy also maintains a Debian package for bedtools
that is likely to be found in its derivatives like Ubuntu. Many thanks to
Charles for doing this.

.. code-block:: bash

    apt-get install bedtools

**Homebrew**. Carlos Borroto has made BEDTools available on the bedtools
package manager for OSX.

.. code-block:: bash

    brew tap homebrew/science
    brew install bedtools

**MacPorts**. Alternatively, the MacPorts ports system can be used to install BEDTools on OSX.

.. code-block:: bash

    port install bedtools

-----------------------------
Development versions
-----------------------------

The development version of bedtools is maintained in a Github
`repository <https://www.github.com/arq5x/bedtools2>`_. Bug fixes are addressed
in this repository prior to release, so there may be situations where you will
want to use a development version of bedtools prior to its being promoted to
a stable release.  One would either clone the repository with ``git``, as
follows and then compile the source code as describe above:

.. code-block:: bash

    git clone https://github.com/arq5x/bedtools2.git


or, one can download the source code as a ``.zip`` file using the Github
website.  Once the zip file is downloaded and uncompressed with the ``unzip``
command, one can compile and install using the instructions above.

    .. image:: images/github-zip-button.png
