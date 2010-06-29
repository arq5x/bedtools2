/*****************************************************************************
  pairToPair.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "pairToPair.h"


/*
	Constructor
*/
PairToPair::PairToPair(string &bedAFilePE, string &bedBFilePE, float &overlapFraction, 
						   string searchType, bool ignoreStrand) {

	_bedAFilePE      = bedAFilePE;
	_bedBFilePE      = bedBFilePE;
	_overlapFraction = overlapFraction;
	_searchType      = searchType;
	_ignoreStrand    = ignoreStrand;
	
	_bedA = new BedFilePE(bedAFilePE);
	_bedB = new BedFilePE(bedBFilePE);
	
	IntersectPairs();
}


/*
	Destructor
*/
PairToPair::~PairToPair(void) {
}



void PairToPair::IntersectPairs() {
	
	// load the "B" bed file into a map so
	// that we can easily compare "A" to it for overlaps
	_bedB->loadBedPEFileIntoMap();
	
	int lineNum = 0;	
	vector<BEDCOV> hitsA1B1, hitsA1B2, hitsA2B1, hitsA2B2;
	// reserve some space
	hitsA1B1.reserve(100); hitsA1B2.reserve(100); hitsA2B1.reserve(100); hitsA2B2.reserve(100);
	
	BedLineStatus bedStatus;
	BEDPE a, nullBedPE;
	
	_bedA->Open();	
	while ((bedStatus = _bedA->GetNextBedPE(a, lineNum)) != BED_INVALID) {
		if (bedStatus == BED_VALID) {
            // identify overlaps b/w the pairs
			FindOverlaps(a, hitsA1B1, hitsA1B2, hitsA2B1, hitsA2B2, _searchType);
		
			// reset space for next BEDPE
			hitsA1B1.clear(); hitsA1B2.clear(); hitsA2B1.clear(); hitsA2B2.clear();		
			a = nullBedPE;
		}
	}
	_bedA->Close();
}
// END IntersectPE



void PairToPair::FindOverlaps(const BEDPE &a, vector<BEDCOV> &hitsA1B1, vector<BEDCOV> &hitsA1B2, 
							  vector<BEDCOV> &hitsA2B1, vector<BEDCOV> &hitsA2B2, string type) {

	// list of hits on each end of BEDPE
	// that exceed the requested overlap fraction
	vector<BEDCOV> qualityHitsA1B1;
	vector<BEDCOV> qualityHitsA1B2;
	vector<BEDCOV> qualityHitsA2B1;
	vector<BEDCOV> qualityHitsA2B2;

	// count of hits on each end of BEDPE
	// that exceed the requested overlap fraction
	int numOverlapsA1B1 = 0;
	int numOverlapsA1B2 = 0;
	int numOverlapsA2B1 = 0;
	int numOverlapsA2B2 = 0;


	// Find the _potential_ hits between each end of A and B
	_bedB->FindOverlapsPerBin(1, a.chrom1, a.start1, a.end1, a.strand1, hitsA1B1, !(_ignoreStrand));	// hits between A1 to B1
	_bedB->FindOverlapsPerBin(1, a.chrom2, a.start2, a.end2, a.strand2, hitsA2B1, !(_ignoreStrand));	// hits between A2 to B1
	_bedB->FindOverlapsPerBin(2, a.chrom1, a.start1, a.end1, a.strand1, hitsA1B2, !(_ignoreStrand));	// hits between A1 to B2
	_bedB->FindOverlapsPerBin(2, a.chrom2, a.start2, a.end2, a.strand2, hitsA2B2, !(_ignoreStrand));	// hits between A2 to B2	


	// Now, reduce to the set of hits on each end of A and B that meet the required overlap fraction and orientation.
	FindQualityHitsBetweenEnds(a, 1, hitsA1B1, qualityHitsA1B1, numOverlapsA1B1);	// quality hits between A1 to B1
	FindQualityHitsBetweenEnds(a, 1, hitsA1B2, qualityHitsA1B2, numOverlapsA1B2);	// quality hits between A2 to B1
	FindQualityHitsBetweenEnds(a, 2, hitsA2B1, qualityHitsA2B1, numOverlapsA2B1);	// quality hits between A1 to B2
	FindQualityHitsBetweenEnds(a, 2, hitsA2B2, qualityHitsA2B2, numOverlapsA2B2);	// quality hits between A2 to B2


	int matchCount1 = 0;	
	int matchCount2 = 0;
	if ((numOverlapsA1B1 > 0) || (numOverlapsA2B2 > 0)) {
		FindHitsOnBothEnds(a, qualityHitsA1B1, qualityHitsA2B2, matchCount1);
	}
	if ((numOverlapsA1B2 > 0) || (numOverlapsA2B1 > 0)) {
		FindHitsOnBothEnds(a, qualityHitsA2B1, qualityHitsA1B2, matchCount2);
	}
	
	
	if ((matchCount1 == 0) && (matchCount2 == 0) && (_searchType == "neither")) {
		_bedA->reportBedPENewLine(a);		
	}
}



void PairToPair::FindQualityHitsBetweenEnds(const BEDPE &a, int end, const vector<BEDCOV> &hits, 
											vector<BEDCOV> &qualityHits, int &numOverlaps) {

	if (end == 1) {
		
		vector<BEDCOV>::const_iterator h = hits.begin();
		vector<BEDCOV>::const_iterator hitsEnd = hits.end();
		for (; h != hitsEnd; ++h) {				
			int s = max(a.start1, h->start);
			int e = min(a.end1, h->end);

			// is there enough overlap (default ~ 1bp)
			if ( ((float)(e-s) / (float)(a.end1 - a.start1)) >= _overlapFraction ) {
				numOverlaps++;
				qualityHits.push_back(*h);
			}
		}
		
	}
	else if (end == 2) {
		
		vector<BEDCOV>::const_iterator h = hits.begin();
		vector<BEDCOV>::const_iterator hitsEnd = hits.end();
		for (; h != hitsEnd; ++h) {				
			int s = max(a.start2, h->start);
			int e = min(a.end2, h->end);
			// is there enough overlap (default ~ 1bp)
			if ( ((float)(e-s) / (float)(a.end2 - a.start2)) >= _overlapFraction ) {
				numOverlaps++;
				qualityHits.push_back(*h);
			}
		}
	}
}


void PairToPair::FindHitsOnBothEnds(const BEDPE &a, const vector<BEDCOV> &qualityHitsEnd1, 
									const vector<BEDCOV> &qualityHitsEnd2, int &matchCount) {
	
	map<unsigned int, vector<BEDCOV>, less<int> > hitsMap;
	
	for (vector<BEDCOV>::const_iterator h = qualityHitsEnd1.begin(); h != qualityHitsEnd1.end(); ++h) {
		hitsMap[h->count].push_back(*h);
		matchCount++;
	}
	for (vector<BEDCOV>::const_iterator h = qualityHitsEnd2.begin(); h != qualityHitsEnd2.end(); ++h) {
		hitsMap[h->count].push_back(*h);
		matchCount++;
	}

	for (map<unsigned int, vector<BEDCOV>, less<unsigned int> >::iterator m = hitsMap.begin(); m != hitsMap.end(); ++m) {
		if (m->second.size() == 2) {
			
			BEDCOV b1 = m->second[0];
			BEDCOV b2 = m->second[1];
			
			if (_searchType == "both") {
				_bedA->reportBedPETab(a);
				printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n", b1.chrom.c_str(), b1.start, b1.end,
																   b2.chrom.c_str(), b2.start, b2.end,
																   b1.name.c_str(), b1.score.c_str(), 
																   b1.strand.c_str(), b2.strand.c_str());
			}
		}
	}
}

