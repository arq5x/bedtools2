/*****************************************************************************
  intersectBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
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

 

void BedIntersect::IntersectBed(istream &bedInput) {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	bedB->loadBedFileIntoMap();

	string bedLine;                                                                                                                    
	int lineNum = 0;					// current input line number
	vector<BED> hits;					// vector of potential hits
	vector<string> bedFields;			// vector for a BED entry
	
	// reserve some space
	hits.reserve(100);
	bedFields.reserve(12);	
		
	// process each entry in A
	while (getline(bedInput, bedLine)) {

		lineNum++;
		Tokenize(bedLine,bedFields);
		BED a;
		
		// find the overlaps with B if it's a valid BED entry. 
		if (bedA->parseLine(a, bedFields, lineNum)) {
			FindOverlaps(a, hits);
			hits.clear();
		}
		
		// reset for the next input line
		bedFields.clear();
	}
}


void BedIntersect::DetermineBedInput() {
	if (bedA->bedFile != "stdin") {   // process a file
		ifstream beds(bedA->bedFile.c_str(), ios::in);
		if ( !beds ) {
			cerr << "Error: The requested bed file (" << bedA->bedFile << ") could not be opened. Exiting!" << endl;
			exit (1);
		}
		IntersectBed(beds);
	}
	else {   						// process stdin
		IntersectBed(cin);		
	}
}