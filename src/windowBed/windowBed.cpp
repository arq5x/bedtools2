// 
//  windowBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Looks for overlaps between features in two BED files.
//

/*
	Includes
*/
#include "windowBed.h"


/*
	Constructor
*/
BedWindow::BedWindow(string &bedAFile, string &bedBFile, int &slop, bool &anyHit, bool &noHit, bool &writeCount) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;

	this->slop = slop;
	this->anyHit = anyHit;
	this->noHit = noHit;
	this->writeCount = writeCount;
	
	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
}

/*
	Destructor
*/
BedWindow::~BedWindow(void) {
}



/*
	reportA
	
	Writes the _original_ BED entry for A.
	Works for BED3 - BED6.
*/
void BedWindow::reportA(const BED &a) {
	
	if (bedA->bedType == 3) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end;
	}
	else if (bedA->bedType == 4) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t"
		<< a.name;
	}
	else if (bedA->bedType == 5) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t"
		<< a.name << "\t" << a.score;
	}
	else if (bedA->bedType == 6) {
		cout << a.chrom << "\t" << a.start << "\t" << a.end << "\t" 
		<< a.name << "\t" << a.score << "\t" << a.strand;
	}
}



/*
	reportB
	
	Writes the _original_ BED entry for B.
	Works for BED3 - BED6.
*/
void BedWindow::reportB(const BED &b) {
	if (bedB->bedType == 3) {
		cout << b.chrom << "\t" << b.start << "\t" << b.end;
	}
	else if (bedB->bedType == 4) {
		cout << b.chrom << "\t" << b.start << "\t" << b.end << "\t"
		<< b.name;
	}
	else if (bedB->bedType == 5) {
		cout << b.chrom << "\t" << b.start << "\t" << b.end << "\t"
		<< b.name << "\t" << b.score;
	}
	else if (bedB->bedType == 6) {
		cout << b.chrom << "\t" << b.start << "\t" << b.end << "\t" 
		<< b.name << "\t" << b.score << "\t" << b.strand;
	}
}



void BedWindow::FindWindowOverlaps(BED &a, vector<BED> &hits) {
	
	// update the current feature's start and end
	// according to the slop requested (slop = 0 by default)
	int aFudgeStart = 0;
	int aFudgeEnd;

	if ((a.start - this->slop) > 0) {
		aFudgeStart = a.start - this->slop;
	}
	else {
		aFudgeStart = 0;
	}
	aFudgeEnd = a.end + this->slop;

	
	bedB->binKeeperFind(bedB->bedMap[a.chrom], aFudgeStart, aFudgeEnd, hits);

	int numOverlaps = 0;
	for (vector<BED>::iterator h = hits.begin(); h != hits.end(); ++h) {
	
		int s = max(aFudgeStart, h->start);
		int e = min(aFudgeEnd, h->end);
	
		if (s < e) {
			// is there enough overlap (default ~ 1bp)
			if ( ((float)(e-s) / (float)(a.end - a.start)) > 0 ) { 
				numOverlaps++;	
				if (!anyHit && !noHit && !writeCount) {			
					reportA(a); cout << "\t";
					reportB(*h); cout << endl;
				}
			}
		}
	}
	if (anyHit && (numOverlaps >= 1)) {
		reportA(a); cout << endl;
	}
	else if (writeCount) {
		reportA(a); cout << "\t" << numOverlaps << endl;
	}
	else if (noHit && (numOverlaps == 0)) {
		reportA(a); cout << endl;
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
			
			vector<string> bedFields;
			Tokenize(bedLine,bedFields);

			lineNum++;
			if (bedA->parseBedLine(a, bedFields, lineNum)) {
				vector<BED> hits;
				FindWindowOverlaps(a, hits);
			}
		}
	}
	// "A" is being passed via STDIN.
	else {
		
		BED a;
		while (getline(cin, bedLine)) {

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
// END WindowIntersectBed


