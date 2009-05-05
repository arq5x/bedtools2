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
						   bool &noHit, bool &writeCount, bool &forceStrand) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;
	this->anyHit = anyHit;
	this->noHit = noHit;
	this->writeA = writeA;	
	this->writeB = writeB;
	this->writeCount = writeCount;
	this->overlapFraction = overlapFraction;
	this->forceStrand = forceStrand;

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
	
	for (vector<BED>::iterator h = hits.begin(); h != hits.end(); ++h) {
	
		// if forcing strandedness, move on if the hit
		// is not on the same strand as A.
		if ((this->forceStrand) && (a.strand != h->strand)) {
			continue;		// continue force the next iteration of the for loop.
		}
	
		int s = max(a.start, h->start);
		int e = min(a.end, h->end);

		if (s < e) {
			
			// is there enough overlap relative to the user's request?
			// (default ~ 1bp)
			if ( ((float)(e-s) / (float)(a.end - a.start)) >= this->overlapFraction ) { 
				
				numOverlaps++;	
				
				if (!anyHit && !noHit && !writeCount) {			
					if (!writeB) {
						if (writeA) {
							bedA->reportBed(a); cout << endl;
						}
						else {
							bedA->reportBedRange(a,s,e);  cout << endl;
						}
					}
					else {
						if (writeA) {
							bedA->reportBed(a); cout << "\t";
							bedB->reportBed(*h); cout << endl;
						}
						else {
							bedA->reportBedRange(a,s,e); cout << "\t";
							bedB->reportBed(*h); cout << endl;										
						}
					}
				}
				
			}
		}
	}
	if (anyHit && (numOverlaps >= 1)) {
		bedA->reportBed(a); cout << endl;
	}
	else if (writeCount) {
		bedA->reportBed(a); cout << "\t" << numOverlaps << endl;
	}
	else if (noHit && (numOverlaps == 0)) {
		bedA->reportBed(a); cout << endl;
	}
}

 

void BedIntersect::IntersectBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();

	string bedLine;
	BED bedEntry;                                                                                                                        
	int lineNum = 0;

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
			
			// split the current line into ditinct fields
			vector<string> bedFields;
			Tokenize(bedLine,bedFields);

			lineNum++;
			
			// find the overlaps with B if it's a valid BED entry. 
			if (bedA->parseBedLine(a, bedFields, lineNum)) {
				vector<BED> hits;
				FindOverlaps(a, hits);
			}
		}
	}
	// "A" is being passed via STDIN.
	else {
		
		BED a;
		// process each entry in A
		while (getline(cin, bedLine)) {

			// split the current line into ditinct fields
			vector<string> bedFields;
			Tokenize(bedLine,bedFields);

			lineNum++;
			
			// find the overlaps with B if it's a valid BED entry. 
			if (bedA->parseBedLine(a, bedFields, lineNum)) {
				vector<BED> hits;
				FindOverlaps(a, hits);
			}
		}
	}
}
// END Intersect


