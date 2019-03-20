#    Copyright (C) 2017 Genome Research Ltd.
#
#    Author: Robert Davies <rmd@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# First field:
#   INIT = initialisation, not counted in testing
#   P = expected to pass (zero return; expected output matches, if present)
#   N = expected to return non-zero
#   F = expected to fail
#
# Second field (P/N/F only):
#   Filename of expected output.  If '.', output is not checked
#
# Rest:
#   Command to execute.  $bgzip and $tabix are replaced with the path to
# bgzip and tabix. 

# TBI index on VCF
INIT $bgzip -c vcf_file.vcf > vcf_file.tbi.tmp.vcf.gz
P . $tabix -f -p vcf vcf_file.tbi.tmp.vcf.gz
P vcf_file.1.3000151.out $tabix vcf_file.tbi.tmp.vcf.gz 1:3000151-3000151
P vcf_file.2.3199812.out $tabix vcf_file.tbi.tmp.vcf.gz 2:3199812-3199812

# CSI index on VCF
INIT $bgzip -c vcf_file.vcf > vcf_file.csi.tmp.vcf.gz
P . $tabix -f -C -p vcf vcf_file.csi.tmp.vcf.gz
P vcf_file.1.3000151.out $tabix vcf_file.csi.tmp.vcf.gz 1:3000151-3000151
P vcf_file.2.3199812.out $tabix vcf_file.csi.tmp.vcf.gz 2:3199812-3199812

# VCF file with chromosome > 2^29-1 bases long
# TBI cannot index this file, so building the index should fail
INIT $bgzip -c large_chr.vcf > large_chr.tmp.vcf.gz
N . $tabix -f -p vcf large_chr.tmp.vcf.gz

# CSI can handle positions > 2^29-1, so building should work
P . $tabix -f -C -p vcf large_chr.tmp.vcf.gz
P large_chr.20.1.2147483647.out $tabix large_chr.tmp.vcf.gz chr20:1-2147483647

# TBI index on BED
INIT $bgzip -c bed_file.bed > bed_file.tbi.tmp.bed.gz
P . $tabix -f -p bed bed_file.tbi.tmp.bed.gz
P bed_file.Y.100200.out $tabix bed_file.tbi.tmp.bed.gz Y:100200-100200

# TBI index on GFF3
INIT $bgzip -c gff_file.gff > gff_file.tbi.tmp.gff.gz
P . $tabix -f -p gff gff_file.tbi.tmp.gff.gz
P gff_file.X.2934832.2935190.out $tabix gff_file.tbi.tmp.gff.gz X:2934832-2935190
