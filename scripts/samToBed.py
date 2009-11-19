#!/usr/bin/env python
# encoding: utf-8
"""
samToBed.py

Created by Aaron Quinlan on 2009-08-27.

PURPOSE:
	Convert SAM format alignment files to BED format, where each alignment
	will conrrespond to a distnct BED line.
"""

import sys
import getopt
import re

help_message = '''

samToBed.py -s <sam> -t <alignment type>

ABSTRACT: Converts aligned reads in SAM format to BED format.
	
OPTIONS:
	-s	The SAM file to be converted to BED (use "stdin" for piped input)
	-t	What types of alignments should be reported?
			"all"	all aligned reads will be reported (Default)
			"con"	only concordant pairs will be reported
			"dis"	only discordant pairs will be reported

EXAMPLE:
	Can be used with samtools to extract alignments and compare them to BED
	annotations.
	
	(1) Land a BED file first.
	$ samtools view reads.sorted.bam > read.sorted.sam
	$ samToBed.py -s reads.sorted.sam -t all > reads.sorted.bed
	$ intersectBed -a reads.sorted.bed -b refseq.bed > reads.intersect.refseq.bed

	(2) "One-liner.
	$ samtools view reads.sorted.bam | samToBed.py -s stdin -t all | \ 
	  intersectBed -a stdin -b refseq.bed > reads.intersect.refseq.bed

'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def processSAM(file, alignType):
	"""
		Read a SAM file (or stdin) and convert each line to BED format.
	"""	
	if (file != "stdin"):	
		for line in open(file,'r'):
			samLine = splitLine(line.strip())
			makeBED(samLine, alignType)
		f.close()
	else:
		for line in sys.stdin:
			samLine = splitLine(line.strip())
			makeBED(samLine, alignType)		
	
					
def makeBED(samFields, aType):
	
	samFlag = int(samFields[1])
	
	# Only create a BED entry if the read was aligned
	if (not (samFlag & 0x0004)):
		
		chrom = samFields[2]
		start = str(int(samFields[3])-1)
		end = str(int(samFields[3]) + len(samFields[9]) - 1)
		name = samFields[0]
		score = samFields[1] 			# we'll use the SAM flag as a score.	
		strand = getStrand(samFlag)

		
		# Write the BED line per user's request.
		printBED(aType, properPairing(samFlag), 
				chrom, start, end, name, score, strand)

		
def splitLine(line, delim="\t"):
	splitline = line.split(delim)
	return splitline		


def properPairing(samFlag):
	return samFlag & (0x0002)


def getStrand(samFlag):
	strand = "+"
	if (samFlag & (0x0010)):	# minus strand if true.
		strand = "-"		
	return strand


def printBED(aType, concordant, chrom, start, end, name, score, strand):
	try:
		if (aType == "all"):
			print chrom + "\t" + start + "\t" + end + "\t" + name + "\t" + score + "\t" + strand
		else:
			if (aType == "con") and concordant:
				print chrom + "\t" + start + "\t" + end + "\t" + name + "\t" + score + "\t" + strand
			elif (aType == "dis") and not concordant:
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
			opts, args = getopt.getopt(argv[1:], "hs:t:", ["help", "sam", "type"])
		except getopt.error, msg:
			raise Usage(help_message)
	
		# option processing
		samFile = ""
		aType = "all"
		for option, value in opts:
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-s", "--sam"):
				samFile = value
			if option in ("-t", "--type"):
				aType = value
	
		if (samFile != "stdin"):
			try:
			   f = open(samFile, 'r')
			except IOError, msg:
				raise Usage(help_message)
					
	except Usage, err:
		print >> sys.stderr, str(err.msg)
		return 2	


	processSAM(samFile, aType)

if __name__ == "__main__":
	sys.exit(main())
