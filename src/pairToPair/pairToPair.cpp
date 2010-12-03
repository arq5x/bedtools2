/*****************************************************************************
  pairToPair.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "pairToPair.h"


/*
	Constructor
*/
PairToPair::PairToPair(string &bedAFilePE, string &bedBFilePE, float &overlapFraction, 
						   string searchType, bool ignoreStrand, int slop, bool strandedSlop) {

	_bedAFilePE      = bedAFilePE;
	_bedBFilePE      = bedBFilePE;
	_overlapFraction = overlapFraction;
	_searchType      = searchType;
	_ignoreStrand    = ignoreStrand;
    _slop            = slop;
    _strandedSlop    = strandedSlop;
	
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
	vector<MATE> hitsA1B1, hitsA1B2, hitsA2B1, hitsA2B2;
	// reserve some space
	hitsA1B1.reserve(100); hitsA1B2.reserve(100); hitsA2B1.reserve(100); hitsA2B2.reserve(100);
	
	BedLineStatus bedStatus;
	BEDPE a, nullBedPE;
	
	_bedA->Open();	
	while ((bedStatus = _bedA->GetNextBedPE(a, lineNum)) != BED_INVALID) {
		if (bedStatus == BED_VALID) {
            // identify overlaps b/w the pairs
			FindOverlaps(a, hitsA1B1, hitsA1B2, hitsA2B1, hitsA2B2);
		
			// reset space for next BEDPE
			hitsA1B1.clear(); hitsA1B2.clear(); hitsA2B1.clear(); hitsA2B2.clear();		
			a = nullBedPE;
		}
	}
	_bedA->Close();
}
// END IntersectPE



void PairToPair::FindOverlaps(const BEDPE &a, 
                              vector<MATE> &hitsA1B1, 
                              vector<MATE> &hitsA1B2, 
							  vector<MATE> &hitsA2B1, 
							  vector<MATE> &hitsA2B2) {

	// list of hits on each end of BEDPE
	// that exceed the requested overlap fraction
	vector<MATE> qualityHitsA1B1;
	vector<MATE> qualityHitsA1B2;
	vector<MATE> qualityHitsA2B1;
	vector<MATE> qualityHitsA2B2;

	// count of hits on each end of BEDPE
	// that exceed the requested overlap fraction
	int numOverlapsA1B1 = 0;
	int numOverlapsA1B2 = 0;
	int numOverlapsA2B1 = 0;
	int numOverlapsA2B2 = 0;
	
	// add the appropriate slop to the starts and ends 
    CHRPOS start1 = a.start1;
    CHRPOS end1   = a.end1;
    CHRPOS start2 = a.start2;
    CHRPOS end2   = a.end2;
    
    if (_strandedSlop == true) {
        if (a.strand1 == "+") 
            end1   += _slop;
        else
            start1 -= _slop;
        if (a.strand2 == "+") 
            end2   += _slop;
        else
            start2 -= _slop;
    }
    else {
        start1 -= _slop;
        start2 -= _slop;
        end1   += _slop;
        end2   += _slop;        
    }

	// Find the _potential_ hits between each end of A and B
	_bedB->FindOverlapsPerBin(1, a.chrom1, start1, end1, a.strand1, hitsA1B1, !(_ignoreStrand));	// hits b/w A1 & B1
	_bedB->FindOverlapsPerBin(1, a.chrom2, start2, end2, a.strand2, hitsA2B1, !(_ignoreStrand));	// hits b/w A2 & B1
	_bedB->FindOverlapsPerBin(2, a.chrom1, start1, end1, a.strand1, hitsA1B2, !(_ignoreStrand));	// hits b/w A1 & B2
	_bedB->FindOverlapsPerBin(2, a.chrom2, start2, end2, a.strand2, hitsA2B2, !(_ignoreStrand));	// hits b/w A2 & B2	

	// Now, reduce to the set of hits on each end of A and B 
	// that meet the required overlap fraction and orientation.
    // FindQualityHitsBetweenEnds(start1, end1, start2, end2, hitsA1B1, qualityHitsA1B1, numOverlapsA1B1);     // quality hits b/w A1 & B1
    // FindQualityHitsBetweenEnds(start1, end1, start2, end2, hitsA1B2, qualityHitsA1B2, numOverlapsA1B2);     // quality hits b/w A1 & B2
    // FindQualityHitsBetweenEnds(start1, end1, start2, end2, hitsA2B1, qualityHitsA2B1, numOverlapsA2B1);     // quality hits b/w A2 & B1
    // FindQualityHitsBetweenEnds(start1, end1, start2, end2, hitsA2B2, qualityHitsA2B2, numOverlapsA2B2);     // quality hits b/w A2 & B2
	FindQualityHitsBetweenEnds(start1, end1, hitsA1B1, qualityHitsA1B1, numOverlapsA1B1);	   // quality hits b/w A1 & B1
	FindQualityHitsBetweenEnds(start1, end1, hitsA1B2, qualityHitsA1B2, numOverlapsA1B2);	   // quality hits b/w A1 & B2
	FindQualityHitsBetweenEnds(start2, end2, hitsA2B1, qualityHitsA2B1, numOverlapsA2B1);	   // quality hits b/w A2 & B1
	FindQualityHitsBetweenEnds(start2, end2, hitsA2B2, qualityHitsA2B2, numOverlapsA2B2);	   // quality hits b/w A2 & B2


	int matchCount1 = 0;	
	int matchCount2 = 0;
	if (_searchType == "neither" || _searchType == "both") {
    	if ((numOverlapsA1B1 > 0) || (numOverlapsA2B2 > 0))
    		FindHitsOnBothEnds(a, qualityHitsA1B1, qualityHitsA2B2, matchCount1);
    	if ((numOverlapsA1B2 > 0) || (numOverlapsA2B1 > 0)) 
    		FindHitsOnBothEnds(a, qualityHitsA2B1, qualityHitsA1B2, matchCount2);
        
        // report the fact that no hits were found iff _searchType is neither.
        if ((matchCount1 == 0) && (matchCount2 == 0) && (_searchType == "neither")) {
    		_bedA->reportBedPENewLine(a);		
    	}
	}
	else if (_searchType == "either") {
	    FindHitsOnEitherEnd(a, qualityHitsA1B1, qualityHitsA2B2, matchCount1);
	    FindHitsOnEitherEnd(a, qualityHitsA2B1, qualityHitsA1B2, matchCount2);
	}
}


void PairToPair::FindQualityHitsBetweenEnds(CHRPOS start, CHRPOS end, const vector<MATE> &hits, 
                                            vector<MATE> &qualityHits, int &numOverlaps) {

	vector<MATE>::const_iterator h       = hits.begin();
	vector<MATE>::const_iterator hitsEnd = hits.end();
	for (; h != hitsEnd; ++h) {				
		int s = max(start, h->bed.start);
		int e = min(end, h->bed.end);

		// is there enough overlap (default ~ 1bp)
		if ( ((float)(e-s) / (float)(end - start)) >= _overlapFraction ) {
			numOverlaps++;
			qualityHits.push_back(*h);
		}
	}
}


void PairToPair::FindHitsOnBothEnds(const BEDPE &a, const vector<MATE> &qualityHitsEnd1, 
									const vector<MATE> &qualityHitsEnd2, int &matchCount) {
	
	map<unsigned int, vector<MATE>, less<int> > hitsMap;
	
	for (vector<MATE>::const_iterator h = qualityHitsEnd1.begin(); h != qualityHitsEnd1.end(); ++h) {
		hitsMap[h->lineNum].push_back(*h);
		matchCount++;
	}
	for (vector<MATE>::const_iterator h = qualityHitsEnd2.begin(); h != qualityHitsEnd2.end(); ++h) {
		hitsMap[h->lineNum].push_back(*h);
		matchCount++;
	}

	for (map<unsigned int, vector<MATE>, less<unsigned int> >::iterator m = hitsMap.begin(); m != hitsMap.end(); ++m) {
		if (m->second.size() == 2) {
			
			MATE b1 = m->second[0];
			MATE b2 = m->second[1];
			
			if (_searchType == "both") {
				_bedA->reportBedPETab(a);
				printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s", b1.bed.chrom.c_str(), b1.bed.start, b1.bed.end,
																   b2.bed.chrom.c_str(), b2.bed.start, b2.bed.end,
																   b1.bed.name.c_str(), b1.bed.score.c_str(), 
																   b1.bed.strand.c_str(), b2.bed.strand.c_str());
				for (size_t i = 0; i < b1.bed.otherFields.size(); ++i)
                    printf("\t%s", b1.bed.otherFields[i].c_str());
                printf("\n");
			}
		}
	}
}


void PairToPair::FindHitsOnEitherEnd(const BEDPE &a, const vector<MATE> &qualityHitsEnd1, 
									const vector<MATE> &qualityHitsEnd2, int &matchCount) {
	
	map<unsigned int, vector<MATE>, less<int> > hitsMap;
	
	for (vector<MATE>::const_iterator h = qualityHitsEnd1.begin(); h != qualityHitsEnd1.end(); ++h) {
		hitsMap[h->lineNum].push_back(*h);
		matchCount++;
	}
	for (vector<MATE>::const_iterator h = qualityHitsEnd2.begin(); h != qualityHitsEnd2.end(); ++h) {
		hitsMap[h->lineNum].push_back(*h);
		matchCount++;
	}

	for (map<unsigned int, vector<MATE>, less<unsigned int> >::iterator m = hitsMap.begin(); m != hitsMap.end(); ++m) {
		if (m->second.size() >= 1) {
			
			if ((m->second.size()) == 2) {
    			MATE b1 = m->second[0];
    			MATE b2 = m->second[1];
			    
    			_bedA->reportBedPETab(a);
    			printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s", b1.bed.chrom.c_str(), b1.bed.start, b1.bed.end,
    															   b2.bed.chrom.c_str(), b2.bed.start, b2.bed.end,
    															   b1.bed.name.c_str(), b1.bed.score.c_str(), 
                                                                   b1.bed.strand.c_str(), b2.bed.strand.c_str());
                for (size_t i = 0; i < b1.bed.otherFields.size(); ++i)
    				printf("\t%s", b1.bed.otherFields[i].c_str());
                printf("\n");
			}
			else {
			    MATE b1 = m->second[0];
			    
    			_bedA->reportBedPETab(a);
    			printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s", b1.bed.chrom.c_str(), b1.bed.start, b1.bed.end,
    															   b1.mate->bed.chrom.c_str(), b1.mate->bed.start, b1.mate->bed.end,
    															   b1.bed.name.c_str(), b1.bed.score.c_str(), 
                                                                   b1.bed.strand.c_str(), b1.mate->bed.strand.c_str());
                for (size_t i = 0; i < b1.bed.otherFields.size(); ++i)
    				printf("\t%s", b1.bed.otherFields[i].c_str());
                printf("\n");		
            }
		}
	}
}
