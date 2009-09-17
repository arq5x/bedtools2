// 
//  windowBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Looks for overlaps between features in two BED files.
//
#include "lineFileUtilities.h"
#include "windowBed.h"


/*
	Constructor
*/
BedWindow::BedWindow(string &bedAFile, string &bedBFile, int &leftSlop, int &rightSlop, bool &anyHit, bool &noHit, bool &writeCount, bool &forceStrand) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;

	this->leftSlop = leftSlop;
	this->rightSlop = rightSlop;

	this->anyHit = anyHit;
	this->noHit = noHit;
	this->writeCount = writeCount;
	this->forceStrand = forceStrand;	
	
	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
}



/*
	Destructor
*/
BedWindow::~BedWindow(void) {
}



void BedWindow::FindWindowOverlaps(BED &a, vector<BED> &hits) {
	
	// update the current feature's start and end
	// according to the slop requested (slop = 0 by default)
	int aFudgeStart = 0;
	int aFudgeEnd;

	if ((a.start - this->leftSlop) > 0) {
		aFudgeStart = a.start - this->leftSlop;
	}
	else {
		aFudgeStart = 0;
	}
	aFudgeEnd = a.end + this->rightSlop;

	
	bedB->binKeeperFind(bedB->bedMap[a.chrom], aFudgeStart, aFudgeEnd, hits);

	int numOverlaps = 0;
	for (vector<BED>::iterator h = hits.begin(); h != hits.end(); ++h) {
	
		// if forcing strandedness, move on if the hit
		// is not on the same strand as A.
		if ((this->forceStrand) && (a.strand != h->strand)) {
			continue;		// continue force the next iteration of the for loop.
		}
	
		int s = max(aFudgeStart, h->start);
		int e = min(aFudgeEnd, h->end);
	
		if (s < e) {
			// is there enough overlap (default ~ 1bp)
			if ( ((float)(e-s) / (float)(a.end - a.start)) > 0 ) { 
				numOverlaps++;	
				if (!anyHit && !noHit && !writeCount) {			
					bedA->reportBed(a); cout << "\t";
					bedB->reportBed(*h); cout << endl;
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

 
void BedWindow::WindowIntersectBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();


	string bedLine;
	BED bedEntry;                                                                                                                        
	int lineNum = 0;

	// open the BED file for reading
	if (bedA->bedFile != "stdin") {

		ifstream bed(bedA->bedFile.c_str(), ios::in);
		if ( !bed ) {
			cerr << "Error: The requested bed file (" <<bedA->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		
		BED a;
		while (getline(bed, bedLine)) {
			
			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {	

				vector<string> bedFields;
				Tokenize(bedLine,bedFields);

				lineNum++;
				if (bedA->parseBedLine(a, bedFields, lineNum)) {
					vector<BED> hits;
					FindWindowOverlaps(a, hits);
				}
			}
		}
	}
	// "A" is being passed via STDIN.
	else {
		
		BED a;
		while (getline(cin, bedLine)) {

			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {	

				vector<string> bedFields;
				Tokenize(bedLine,bedFields);

				lineNum++;
				if (bedA->parseBedLine(a, bedFields, lineNum)) {
					vector<BED> hits;
					FindWindowOverlaps(a, hits);
				}
			}
		}
	}
}
// END WindowIntersectBed


