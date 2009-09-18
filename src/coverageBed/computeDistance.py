#!/usr/bin/env python

import sys
import intervals

"""
chr1	65815194	65818190	contig5	2995	chr1	65824190	65826678
chr1	173492516	173494363	contig14	1846	chr1	173471151	173491471
chr1	173470227	173479781	contig22	9553	chr1	173471151	173491471
"""

for line in sys.stdin:
	l = line.strip().split("\t")
	
	s1 = int(l[1])
	e1 = int(l[2])
	s2 = int(l[6])
	e2 = int(l[7])
	
	print l[0] + "\t" + str(s1) + "\t" + str(e1) + "\t" + l[3] + "\t" + str(l[4]) + "\t" + str(intervals.overlap(s1,e1,s2,e2))
