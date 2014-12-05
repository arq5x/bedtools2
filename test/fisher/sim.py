import sys
from random import randint
from subprocess import check_output

nA = 3000
minA, maxA = (20, 5250)
#minA, maxA = (200, 250)

bIntervals = [(x[0], int(x[1]), int(x[2])) for x in (l.split("\t") for l in
    open('hg19.knownCanonical.bed')) if x[0] == "chr1" ]
bIntervals.sort()
genome_size = max(b[2] for b in bIntervals) + 50000

with open('tbb.bed', 'w') as fh:
    for chrom, start, end in bIntervals:
        fh.write("\t".join((chrom, str(start), str(end))) + "\n")

with open('taa.bed', 'w') as fh:
    vals = []
    for i in range(nA):
        s = randint(0, genome_size - maxA)
        e = randint(s + minA, min(genome_size, s + maxA))
        vals.append((s, e))
    for s, e in sorted(vals):
        fh.write("chr1\t%i\t%i\n" % (s, e))
    fh.flush()

print >> open('tgg.genome', 'w'), ("chr1\t%i" % genome_size)

# NOTE: add -m here to make merged output
print check_output("../../bin/bedtools fisher -a taa.bed -b tbb.bed -g tgg.genome", shell=True).strip()
