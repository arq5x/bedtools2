import sys

try:
    xrange
except:
    xrange = range

size = 5000000000
linelen = 200
extra = "ACTGACCCCGAGACGTTTGCATCCTGCACAGCTAGAGATCCTTTATTAAAAGCACACTGT"

with open("bigx.fasta", "wb") as fh:
    fh.write(">chrbig\n")
    line = ("N" * linelen) + "\n"
    for i in xrange(0, size, len(line) - 1):
        fh.write(line)
    fh.write(extra + "\n")
fh.close()

with open("bigx.fasta.fai", "wb") as fh:
    fh.write("chrbig\t%d\t8\t%d\t%d\n" % (size + len(extra), linelen, linelen+1))

