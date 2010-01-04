/*****************************************************************************
  pairToPair.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#ifndef PAIRTOPAIR_H
#define PAIRTOPAIR_H

#include "bedFile.h"
#include "bedFilePE.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class PairToPair {

public:

	// constructor 
	PairToPair(string &, string &, float &, string, bool);

	// destructor
	~PairToPair(void);

 	void IntersectPairs(istream &bedInput);

	void FindOverlaps(const BEDPE &a, vector<BED> &hitsA1B1, vector<BED> &hitsA1B2, 
		vector<BED> &hitsA2B1, vector<BED> &hitsA2B2, string type);

	void FindQualityHitsBetweenEnds(const BEDPE &a, int end, const vector<BED> &hits, 
		vector<BED> &qualityHits, int &numOverlaps);
	
	void FindHitsOnBothEnds(const BEDPE &a, const vector<BED> &qualityHitsEnd1, 
		const vector<BED> &qualityHitsEnd2, int &matchCount);
	
	void DetermineBedPEInput();
		
private:

	string bedAFilePE;
	string bedBFilePE;
	
	float overlapFraction;
	string searchType;
	bool ignoreStrand;

	// instance of a paired-end bed file class.
	BedFilePE *bedA;

	// instance of a bed file class.
	BedFilePE *bedB;
};

#endif /* PAIRTOPAIR_H */
