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

samToBed -s <sam> -t <alignment type>
	
OPTIONS:
	-s	The SAM file to be converted to BED
	-t	What types of alignments should be reported?
			"all"	all aligned reads will be reported (Default)
			"con"	only concordant pairs will be reported
			"dis"	only discordant pairs will be reported

'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def processSAM(file, alignType):
	"""
		Load a SAM file and convert each line to BED format.
		
		We avoid readlines() in this case, as SAM files can 
		be HUGE, and thus loading it into memory could be painful.
	"""		
	for line in open(file,'r'):
		samLine = splitLine(line.strip())
		makeBED(samLine, alignType)
	f.close()	
	
					
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
				
		try:
		   f = open(samFile, 'r')
		except IOError, msg:
			raise Usage(help_message)
				
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		return 2	

	# make a BED file of the SAM file.
	processSAM(samFile, aType)

if __name__ == "__main__":
	sys.exit(main())
