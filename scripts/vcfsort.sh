#!/bin/bash

[ $# -eq 0 ] && { echo "Sorts a VCF file in natural chromosome order";\
                  echo "Usage: $0 [my.vcf | my.vcf.gz]"; exit 1;
                 }
if LC_ALL=C (zless $1 | grep ^#; zless $1 | grep -v ^# | sort -k1,1V -k2,2n);
then
    exit 0
else
    printf 'sort failed. Does your version of sort support the -V option?\n'
    printf 'If not, you should update sort with the latest from GNU coreutils:\n'
    printf 'git clone git://git.sv.gnu.org/coreutils'
fi
