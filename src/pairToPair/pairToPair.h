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

 	void IntersectPairs();

	void FindOverlaps(const BEDPE &a, vector<BEDCOV> &hitsA1B1, vector<BEDCOV> &hitsA1B2, 
		vector<BEDCOV> &hitsA2B1, vector<BEDCOV> &hitsA2B2, string type);

	void FindQualityHitsBetweenEnds(const BEDPE &a, int end, const vector<BEDCOV> &hits, 
		vector<BEDCOV> &qualityHits, int &numOverlaps);
	
	void FindHitsOnBothEnds(const BEDPE &a, const vector<BEDCOV> &qualityHitsEnd1, 
		const vector<BEDCOV> &qualityHitsEnd2, int &matchCount);
	
		
private:

	string _bedAFilePE;
	string _bedBFilePE;
	
	float _overlapFraction;
	string _searchType;
	bool _ignoreStrand;

	// instance of a paired-end bed file class.
	BedFilePE *_bedA;

	// instance of a bed file class.
	BedFilePE *_bedB;
};

#endif /* PAIRTOPAIR_H */
