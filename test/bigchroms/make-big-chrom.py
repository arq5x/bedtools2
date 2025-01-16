import sys

try:
    xrange
except:
    xrange = range

size = 5000000000
linelen = 200
extra = "ACTGACCCCGAGACGTTTGCATCCTGCACAGCTAGAGATCCTTTATTAAAAGCACACTGT"

with open("bigx.fasta", "wb") as fh:
    fh.write(b">chrbig\n")
    line = ("N" * linelen) + "\n"
    for i in xrange(0, size, len(line) - 1):
        fh.write(bytes(line, "utf8"))
    fh.write(bytes(extra + "\n", "utf8"))
fh.close()

with open("bigx.fasta.fai", "wb") as fh:
    fh.write(b"chrbig\t%d\t8\t%d\t%d\n" % (size + len(extra), linelen, linelen+1))

