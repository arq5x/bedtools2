// 
//  subtractBed.cpp
//  BEDTools
//  
//  Created by Aaron Quinlan Spring 2009.
//  Copyright 2009 Aaron Quinlan. All rights reserved.
//
//  Summary:  Removes overlapping segments from a BED entry.
//
#include "lineFileUtilities.h"
#include "subtractBed.h"


/*
	Constructor
*/
BedSubtract::BedSubtract(string &bedAFile, string &bedBFile, float &overlapFraction, bool &forceStrand) {

	this->bedAFile = bedAFile;
	this->bedBFile = bedBFile;
	this->overlapFraction = overlapFraction;
	this->forceStrand = forceStrand;

	this->bedA = new BedFile(bedAFile);
	this->bedB = new BedFile(bedBFile);
}


/*
	Destructor
*/
BedSubtract::~BedSubtract(void) {
}



void BedSubtract::FindOverlaps(BED &a, vector<BED> &hits) {
	
	// find all of the overlaps between a and B.
	bedB->binKeeperFind(bedB->bedMap[a.chrom], a.start, a.end, hits);
	
	//  is A completely spanned by an entry in B?
	//  if so, A should not be reported.
	int numConsumedByB = 0; 
	
	for (vector<BED>::iterator h = hits.begin(); h != hits.end(); ++h) {
		
		// if forcing strandedness, move on if the hit
		// is not on the same strand as A.
		if ((this->forceStrand) && (a.strand != h->strand)) {
			continue;		// continue force the next iteration of the for loop.
		}
	
		int s = max(a.start, h->start);
		int e = min(a.end, h->end);

		if (s < e) {
			
			// is there enough overlap (default ~ 1bp)
			float overlap = ((float)(e-s) / (float)(a.end - a.start));
			if (overlap >= 1.0) {
				numConsumedByB++;
			}
		}
	}
	
	// report the subtraction with each entry in B
	// if A was not "consumed" by any entry in B
	if (numConsumedByB == 0) {
		
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
			
				numOverlaps++;
			
				// is there enough overlap (default ~ 1bp)
				float overlap = ((float)(e-s) / (float)(a.end - a.start));
			
				if ( overlap >= this->overlapFraction ) { 
				
					if (overlap < 1.0) {			// if overlap = 1, then B consumes A and therefore,
													// we won't report A.
						// A	++++++++++++
						// B        ----
						// Res. ====    ====					
						if ( (h->start > a.start) && (h->end < a.end) ) {
							bedA->reportBedRange(a,a.start,h->start); cout << "\n";
							bedA->reportBedRange(a,h->end,a.end); cout << "\n";
						}
						// A	++++++++++++
						// B    ----------
						// Res.           ==        			
						else if (h->start == a.start) {
							bedA->reportBedRange(a,h->end,a.end); cout << "\n";
						}	
						// A	      ++++++++++++
						// B    ----------
						// Res.       ====
						else if (h->start < a.start) {
							bedA->reportBedRange(a,h->end,a.end); cout << "\n";					
						}
						// A	++++++++++++
						// B           ----------
						// Res. =======
						else if (h->start > a.start) {
							bedA->reportBedRange(a,a.start,h->start); cout << "\n";					
						}
					}
				}
			}
		
		}
		if (numOverlaps == 0) {
			// no overlap found, so just report A as-is.
			bedA->reportBed(a); cout << "\n";
		}
	}
}

 

void BedSubtract::SubtractBed() {

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
			
			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {	
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
	// "A" is being passed via STDIN.
	else {
		
		BED a;
		// process each entry in A
		while (getline(cin, bedLine)) {

			if ((bedLine.find("track") != string::npos) || (bedLine.find("browser") != string::npos)) {
				continue;
			}
			else {	
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
}
// END Intersect


