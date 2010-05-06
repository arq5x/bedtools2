/*****************************************************************************
  closestBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "closestBed.h"

const int MAXSLOP = 256000000;  // 2*MAXSLOP = 512 megabases.
			        			// We don't want to keep looking if we
			        			// can't find a nearby feature within 512 Mb.
const int SLOPGROWTH = 2048000;


/*
	Constructor
*/
BedClosest::BedClosest(string &bedAFile, string &bedBFile, bool &forceStrand, string &tieMode) {

	_bedAFile = bedAFile;
	_bedBFile = bedBFile;
	_forceStrand = forceStrand;
	_tieMode = tieMode;

	_bedA = new BedFile(bedAFile);
	_bedB = new BedFile(bedBFile);
	
	FindClosestBed();
}

/*
	Destructor
*/
BedClosest::~BedClosest(void) {
}


/*
	reportNullB
	
	Writes a NULL B entry for cases where no closest BED was found
	Works for BED3 - BED6.
*/
void BedClosest::reportNullB() {
	if (_bedB->bedType == 3) {
		printf("none\t-1\t-1\n");
	}
	else if (_bedB->bedType == 4) {
		printf("none\t-1\t-1\t-1\n");
	}
	else if (_bedB->bedType == 5) {
		printf("none\t-1\t-1\t-1\t-1\n");
	}
	else if (_bedB->bedType == 6) {
		printf("none\t-1\t-1\t-1\t-1\t-1\n");
	}
}




void BedClosest::FindWindowOverlaps(BED &a, vector<BED> &hits) {
	
	int slop = 0;  // start out just looking for overlaps 
		       	   // within the current bin (~128Kb)	

	// update the current feature's start and end

	int aFudgeStart = 0;
	int aFudgeEnd;
	int numOverlaps = 0;
	vector<BED> closestB;
	float maxOverlap = 0;
	int minDistance = 999999999;


	if(_bedB->bedMap.find(a.chrom) != _bedB->bedMap.end()) {

		while ((numOverlaps == 0) && (slop <= MAXSLOP)) {
		
			// add some slop (starting at 0 bases) to a in hopes
			// of finding a hit in B
			if ((a.start - slop) > 0) aFudgeStart = a.start - slop;
			else aFudgeStart = 0;
			
			if ((a.start + slop) < 2 * MAXSLOP) aFudgeEnd = a.end + slop;
			else aFudgeEnd = 2 * MAXSLOP;
		
			_bedB->FindOverlapsPerBin(a.chrom, aFudgeStart, aFudgeEnd, a.strand, hits, _forceStrand);
	
			vector<BED>::const_iterator h = hits.begin();
			vector<BED>::const_iterator hitsEnd = hits.end();
			for (; h != hitsEnd; ++h) {
						
				numOverlaps++;

				// do the actual features overlap?		
				int s = max(a.start, h->start);
				int e = min(a.end, h->end);
				int overlapBases = (e - s);				// the number of overlapping bases b/w a and b
				int aLength = (a.end - a.start);		// the length of a in b.p.
				
				if (s < e) {
					// is there enough overlap (default ~ 1bp)
					float overlap = (float) overlapBases / (float) aLength;
					if ( overlap > 0 ) {			
						// is this hit the closest?
						if (overlap > maxOverlap) {
							closestB.clear();
							closestB.push_back(*h);
							maxOverlap = overlap;
						}
						else if (overlap == maxOverlap) closestB.push_back(*h);
					}
				}
				else if (h->end < a.start){
					if ((a.start - h->end) < minDistance) {
						closestB.clear();
						closestB.push_back(*h);
						minDistance = a.start - h->end;
					}
					else if ((a.start - h->end) == minDistance) closestB.push_back(*h);
				}
				else {
					if ((h->start - a.end) < minDistance) {
						closestB.clear();
						closestB.push_back(*h);
						minDistance = h->start - a.end;
					}
					else if ((h->start - a.end) == minDistance) closestB.push_back(*h);	
				}
			}
			// if no overlaps were found, we'll widen the range 
			// by SLOPGROWTH in each direction and search again.
			slop += SLOPGROWTH;
		}
	}
	else {
		_bedA->reportBedTab(a);
		reportNullB(); 
	}

	if (numOverlaps > 0) {
		
		if (closestB.size() == 1) {		
			_bedA->reportBedTab(a); 
			_bedB->reportBedNewLine(closestB[0]);
		}
		else {
			if (_tieMode == "all") {
				for (vector<BED>::iterator b = closestB.begin(); b != closestB.end(); ++b) {
					_bedA->reportBedTab(a); 
					_bedB->reportBedNewLine(*b);
				}
			}
			else if (_tieMode == "first") {
				_bedA->reportBedTab(a); 
				_bedB->reportBedNewLine(closestB[0]);
			}
			else if (_tieMode == "last") {
				_bedA->reportBedTab(a); 
				_bedB->reportBedNewLine(closestB[closestB.size()-1]);
			}
		}
	}
}

 
void BedClosest::FindClosestBed() {

	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	_bedB->loadBedFileIntoMap();
	
	BED a, nullBed;
	int lineNum = 0;					// current input line number
	vector<BED> hits;					// vector of potential hits
	hits.reserve(100);
	BedLineStatus bedStatus;

	_bedA->Open();
	// process each entry in A in search of the closest feature in B
	bedStatus = _bedA->GetNextBed(a, lineNum);
	while (bedStatus != BED_INVALID) {
		if (bedStatus == BED_VALID) {
			FindWindowOverlaps(a, hits);
			hits.clear();
			a = nullBed;
		}
		bedStatus = _bedA->GetNextBed(a, lineNum);
	}
	_bedA->Close();
}
// END ClosestBed

