############
Installation
############

BEDTools is intended to run in a "command line" environment on UNIX, LINUX and Apple OS X
operating systems. Installing BEDTools involves downloading the latest source code archive followed by
compiling the source code into binaries on your local system. The following commands will install
BEDTools in a local directory on a NIX or OS X machine. Note that the **"<version>"** refers to the
latest posted version number on http://bedtools.googlecode.com/.

Note: *The BEDTools "makefiles" use the GCC compiler. One should edit the Makefiles accordingly if
one wants to use a different compiler.*::

  curl http://bedtools.googlecode.com/files/BEDTools.<version>.tar.gz > BEDTools.tar.gz
  tar -zxvf BEDTools.tar.gz
  cd BEDTools-<version>
  make clean
  make all
  ls bin
  
At this point, one should copy the binaries in BEDTools/bin/ to either usr/local/bin/ or some
other repository for commonly used UNIX tools in your environment. You will typically require
administrator (e.g. "root" or "sudo") privileges to copy to usr/local/bin/. If in doubt, contact you
system administrator for help.

