// 
//  peIntersectBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Looks for overlaps between paired-end reads / SVs (BEDPE format) and a BED file.
//
#include "lineFileUtilities.h"
#include "peIntersectBed.h"


/*
	Constructor
*/

BedIntersectPE::BedIntersectPE(string &bedAFilePE, string &bedBFile, float &overlapFraction, 
						   string &searchType) {

	this->bedAFilePE = bedAFilePE;
	this->bedBFile = bedBFile;
	this->overlapFraction = overlapFraction;
	this->searchType = searchType;
	
	this->bedA = new BedFilePE(bedAFilePE);
	this->bedB = new BedFile(bedBFile);
}


/*
	Destructor
*/

BedIntersectPE::~BedIntersectPE(void) {
}



void BedIntersectPE::FindOverlaps(BEDPE &a, vector<BED> &hits1, vector<BED> &hits2, string &type) {

	// list of hits on each end of BEDPE
	// that exceed the requested overlap fraction
	vector<BED> qualityHits1;
	vector<BED> qualityHits2;

	// count of hits on each end of BEDPE
	// that exceed the requested overlap fraction
	int numOverlapsEnd1 = 0;
	int numOverlapsEnd2 = 0;


	/* 
	Find the quality hits between ***end1*** of the BEDPE and the B BED file
	*/
	bedB->binKeeperFind(bedB->bedMap[a.chrom1], a.start1, a.end1, hits1);
	
	for (vector<BED>::iterator h = hits1.begin(); h != hits1.end(); ++h) {
	
		int s = max(a.start1, h->start);
		int e = min(a.end1, h->end);

		if (s < e) {
			
			// is there enough overlap (default ~ 1bp)
			if ( ((float)(e-s) / (float)(a.end1 - a.start1)) >= this->overlapFraction ) { 
				numOverlapsEnd1++;
				
				if (type == "either") {
					bedA->reportBedPE(a); cout << "\t";
					bedB->reportBed(*h); cout << endl;
				}
				else {
					qualityHits1.push_back(*h);
				}	
			}
		}
	}
	
	
	/* 
	Now find the quality hits between ***end2*** of the BEDPE and the B BED file
	*/
	bedB->binKeeperFind(bedB->bedMap[a.chrom2], a.start2, a.end2, hits2);
	
	for (vector<BED>::iterator h = hits2.begin(); h != hits2.end(); ++h) {
	
		int s = max(a.start2, h->start);
		int e = min(a.end2, h->end);

		if (s < e) {
			
			// is there enough overlap (default ~ 1bp)
			if ( ((float)(e-s) / (float)(a.end2 - a.start2)) >= this->overlapFraction ) { 
				numOverlapsEnd2++;
				
				if (type == "either") {
					bedA->reportBedPE(a); cout << "\t";
					bedB->reportBed(*h); cout << endl;
				}
				else {
					qualityHits2.push_back(*h);
				}
			}
		}
	}
	
	
	/*
		Now report the hits depending on what the user has requested.
	*/
	if (type == "neither") {
		if ( (numOverlapsEnd1 == 0) && (numOverlapsEnd2 == 0) ) {
			bedA->reportBedPE(a); cout << endl;
		}
	}
	else if (type == "xor") {
		if ( (numOverlapsEnd1 > 0) && (numOverlapsEnd2 == 0) ) {
			for (vector<BED>::iterator q = qualityHits1.begin(); q != qualityHits1.end(); ++q) {
				bedA->reportBedPE(a); cout << "\t";
				bedB->reportBed(*q); cout << endl;
			}
		}
		else if ( (numOverlapsEnd1 == 0) && (numOverlapsEnd2 > 0) ) {
			for (vector<BED>::iterator q = qualityHits2.begin(); q != qualityHits2.end(); ++q) {
				bedA->reportBedPE(a); cout << "\t";
				bedB->reportBed(*q); cout << endl;
			}
		}
	}
	else if (type == "both") {
		if ( (numOverlapsEnd1 > 0) && (numOverlapsEnd2 > 0) ) {
			for (vector<BED>::iterator q = qualityHits1.begin(); q != qualityHits1.end(); ++q) {
				bedA->reportBedPE(a); cout << "\t";
				bedB->reportBed(*q); cout << endl;
			}
			for (vector<BED>::iterator q = qualityHits2.begin(); q != qualityHits2.end(); ++q) {
				bedA->reportBedPE(a); cout << "\t";
				bedB->reportBed(*q); cout << endl;
			}
		}
	}
}





void BedIntersectPE::FindSpanningOverlaps(BEDPE &a, vector<BED> &hits, string &type) {

	// count of hits on _between_ end of BEDPE
	// that exceed the requested overlap fraction
	int numOverlaps = 0;


	/* 
	Find the hits between end1 and start2 of the BEDPE and the B BED file
	
	In other words, find the hits between the "span" of the pair
	*/
	
	int spanStart = 0;
	int spanEnd = 0;
	int spanLength = 0;
	if (type == "inspan") {
		spanStart = a.end1;
		spanEnd = a.start2;
	}
	else if (type == "outspan") {
		spanStart = a.start1;
		spanEnd = a.end2;		
	}
	spanLength = spanEnd - spanStart;
	
	bedB->binKeeperFind(bedB->bedMap[a.chrom1], spanStart, spanEnd, hits);
	
	for (vector<BED>::iterator h = hits.begin(); h != hits.end(); ++h) {
	
		int s = max(spanStart, h->start);
		int e = min(spanEnd, h->end);

		// overlap if s < e
		if (s < e) {
			
			// is there enough overlap (default ~ 1bp)
			if ( (float)(e-s) / (float) spanLength >= this->overlapFraction ) { 
				numOverlaps++;

				bedA->reportBedPE(a); cout << "\t";
				bedB->reportBed(*h); cout << endl;
			}
		}
	}
}


 

void BedIntersectPE::IntersectBedPE() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();

	string bedLine;                                                                                                                     
	int lineNum = 0;

	// are we dealing with a file?
	if (bedA->bedFile != "stdin") {

		// open the BEDPE file for reading
		ifstream bed(bedA->bedFile.c_str(), ios::in);
		if ( !bed ) {
			cerr << "Error: The requested bed file (" <<bedA->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		
		BEDPE a;
		// process each entry in A
		while (getline(bed, bedLine)) {
			
			// split the current line into ditinct fields
			vector<string> bedFields;
			Tokenize(bedLine,bedFields);

			lineNum++;
			
			// find the overlaps with B if it's a valid BED entry. 
			if (bedA->parseBedPELine(a, bedFields, lineNum)) {
				
				if ((this->searchType == "inspan") || (this->searchType == "outspan")) {
					vector<BED> hits;
					if (a.chrom1 == a.chrom2) {
						FindSpanningOverlaps(a, hits, this->searchType);
					}
				}
				else {
					vector<BED> hits1, hits2;
					FindOverlaps(a, hits1, hits2, this->searchType);
				}
				
			}
		}
	}
	// "A" is being passed via STDIN.
	else {
		
		BEDPE a;
		// process each entry in A
		while (getline(cin, bedLine)) {

			// split the current line into ditinct fields
			vector<string> bedFields;
			Tokenize(bedLine,bedFields);

			lineNum++;
			
			// find the overlaps with B if it's a valid BED entry. 
			if (bedA->parseBedPELine(a, bedFields, lineNum)) {

				if ((this->searchType == "inspan") || (this->searchType == "outspan")) {
					vector<BED> hits;
					if (a.chrom1 == a.chrom2) {
						FindSpanningOverlaps(a, hits, this->searchType);
					}
				}
				else {
					vector<BED> hits1, hits2;
					FindOverlaps(a, hits1, hits2, this->searchType);
				}
			}
		}
	}
}
// END IntersectPE


