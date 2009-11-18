#!/usr/bin/env python
# encoding: utf-8
"""
gffToBed.py

Created by Aaron Quinlan on 2009-08-27.

PURPOSE:
	Convert GFF files to BED format.
"""

import sys
import getopt
import re

help_message = '''

gffToBed -g <gff>
	
OPTIONS:
'''

"""
GFF Definition:

## Header lines

0. (used) seqname - The name of the sequence. Must be a chromosome or scaffold.
1. (ignd) source - The program that generated this feature.
2. (used) feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
3. (used) start - The starting position of the feature in the sequence. The first base is numbered 1.
4. (used) end - The ending position of the feature (inclusive).
5. (used) score - A score between 0 and 1000. 
6. (used) strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
7. (ignd) frame - If the feature is a coding exon, 
8. (ignd) group - All lines with the same group are linked together into a single item.
"""

class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def processGFF(file):
	"""
		Load a GFF file and convert each line to BED format.
	"""		
	for line in open(file,'r'):
		gffLine = splitLine(line.strip())
		
		if (gffLine[:2] != "##"):
			makeBED(gffLine)
	f.close()	
	
					
def makeBED(gffFields):
	"""
		Convert a single GFF line to BED.
	"""			

	chrom = gffFields[0]
	start = str(int(gffFields[3])-1)
	end = str(int(gffFields[4]))
	name = str(gffFields[2])
	score = str(getScore(gffFields[5])) 		# Default to 0 if "."
	strand = str(getStrand(gffFields[6]))	# Default to + if "."

	# Write the BED line per user's request.
	printBED(chrom, start, end, name, score, strand)

		
def splitLine(line, delim="\t"):
	splitline = line.split(delim)
	return splitline		


def getScore(gffScore):
	if (gffScore == "."):
		return 0
	else:
		return gffScore


def getStrand(gffStrand):
	if (gffStrand == "."):
		return "+"
	else:
		return gffStrand	


def printBED(chrom, start, end, name, score, strand):
	try:
		print chrom + "\t" + start + "\t" + end + "\t" + name + "\t" + score + "\t" + strand
	except IOError, e:
		sys.exit()
	except KeyboardInterrupt, e:
		sys.exit()	
		
		
def main(argv=None):
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hg:", ["help", "gff"])
		except getopt.error, msg:
			raise Usage(help_message)
	
		# option processing
		gffFile = ""
		for option, value in opts:
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-g", "--gff"):
				gffFile = value
				
		try:
		   f = open(gffFile, 'r')
		except IOError, msg:
			raise Usage(help_message)
				
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		return 2	

	# make a BED file of the GFF file.
	processGFF(gffFile)

if __name__ == "__main__":
	sys.exit(main())
