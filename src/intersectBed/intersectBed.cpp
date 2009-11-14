// 
//  intersectBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Looks for overlaps between features in two BED files.
//
#include "lineFileUtilities.h"
#include "intersectBed.h"


/*
	Constructor
*/
BedIntersect::BedIntersect(string &bedAFile, string &bedBFile, bool &anyHit, 
						   bool &writeA, bool &writeB, float &overlapFraction, 
						   bool &noHit, bool &writeCount, bool &forceStrand, bool &reciprocal) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;
	this->anyHit = anyHit;
	this->noHit = noHit;
	this->writeA = writeA;	
	this->writeB = writeB;
	this->writeCount = writeCount;
	this->overlapFraction = overlapFraction;
	this->forceStrand = forceStrand;
	this->reciprocal = reciprocal;

	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
}


/*
	Destructor
*/
BedIntersect::~BedIntersect(void) {
}


void BedIntersect::FindOverlaps(BED &a, vector<BED> &hits) {
	
	// find all of the overlaps between a and B.
	bedB->binKeeperFind(bedB->bedMap[a.chrom], a.start, a.end, hits);

	int numOverlaps = 0;
	
	// should we print each overlap, or does the 
	// user want summary information?
	bool printable = true;			
	if (anyHit || noHit || writeCount) {
		printable = false;
	}
	
	for (vector<BED>::const_iterator h = hits.begin(); h != hits.end(); ++h) {
	
		// if forcing strandedness, move on if the hit
		// is not on the same strand as A.
		if ((this->forceStrand) && (a.strand != h->strand)) {
			continue;		// continue force the next iteration of the for loop.
		}
		else {
			int s = max(a.start, h->start);
			int e = min(a.end, h->end);
			
			// is there enough overlap relative to the user's request?
			// (default ~ 1bp)
			if ( ((float)(e-s) / (float)(a.end - a.start)) >= this->overlapFraction ) { 
			
				// Report the hit if the user doesn't care about reciprocal overlap
				// between A and B.
				if (!reciprocal) {
				
					numOverlaps++;		// we have another hit for A
				
					if (!writeB && printable) {
						if (writeA) {
							bedA->reportBedNewLine(a);
						}
						else {
							bedA->reportBedRangeNewLine(a,s,e);
						}
					}
					else if (printable) {
						if (writeA) {
							bedA->reportBedTab(a);
							bedB->reportBedNewLine(*h);
						}
						else {
							bedA->reportBedRangeTab(a,s,e);
							bedB->reportBedNewLine(*h);									
						}
					}
				}
				else {
				
					float bOverlap = ((float)(e-s) / (float)(h->end - h->start));
				
					if (bOverlap >= this->overlapFraction) {
					
						numOverlaps++;		// we have another hit for A
					
						if (!writeB && printable) {
							if (writeA) {
								bedA->reportBedNewLine(a);
							}
							else {
								bedA->reportBedRangeNewLine(a,s,e);
							}
						}
						else if (printable) {
							if (writeA) {
								bedA->reportBedTab(a);
								bedB->reportBedNewLine(*h);
							}
							else {
								bedA->reportBedRangeTab(a,s,e);
								bedB->reportBedNewLine(*h);									
							}
						}
					}
				}
			}
		}
	}
	if (anyHit && (numOverlaps >= 1)) {
		bedA->reportBedNewLine(a);
	}
	else if (writeCount) {
		bedA->reportBedTab(a); 
		printf("%d\n", numOverlaps);
	}
	else if (noHit && (numOverlaps == 0)) {
		bedA->reportBedNewLine(a);
	}
}

 

void BedIntersect::IntersectBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();

	string bedLine;
	BED bedEntry;                                                                                                                        
	int lineNum = 0;
	vector<BED> hits;					// vector or potential hits
	vector<string> bedFields;	// vector for a BED entry
	
	// reserve some space
	hits.reserve(100);
	bedFields.reserve(6);	
	
	// are we dealing with a file?
	if (bedA->bedFile != "stdin") {

		// open the BED file for reading
		ifstream bed(bedA->bedFile.c_str(), ios::in);
		if ( !bed ) {
			cerr << "Error: The requested bed file (" <<bedA->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		
		BED a;
		// process each entry in A
		while (getline(bed, bedLine)) {
	
			lineNum++;
			
			if (lineNum > 1) {
				Tokenize(bedLine,bedFields);

				// find the overlaps with B if it's a valid BED entry. 
				if (bedA->parseBedLine(a, bedFields, lineNum)) {
					FindOverlaps(a, hits);
					hits.clear();
				}
				bedFields.clear();
			}
			else {
				if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
					continue;
				}
				else {
					Tokenize(bedLine,bedFields);
					// find the overlaps with B if it's a valid BED entry. 
					if (bedA->parseBedLine(a, bedFields, lineNum)) {
						FindOverlaps(a, hits);
						hits.clear();
					}
					bedFields.clear();
				}
			}
		}
	}
	// "A" is being passed via STDIN.
	else {
		
		BED a;
		// process each entry in A
		while (getline(cin, bedLine)) {

			lineNum++;
			
			if (lineNum > 1) {
				Tokenize(bedLine,bedFields);

				// find the overlaps with B if it's a valid BED entry. 
				if (bedA->parseBedLine(a, bedFields, lineNum)) {
					FindOverlaps(a, hits);
					hits.clear();
				}
				bedFields.clear();
			}
			else {
				if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
					continue;
				}
				else {
					Tokenize(bedLine,bedFields);
					// find the overlaps with B if it's a valid BED entry. 
					if (bedA->parseBedLine(a, bedFields, lineNum)) {
						FindOverlaps(a, hits);
						hits.clear();
					}
					bedFields.clear();
				}
			}
		}
	}
}
// END Intersect


