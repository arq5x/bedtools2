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
	PairToPair(string &, string &, float &, string &);

	// destructor
	~PairToPair(void);

 	void IntersectPairs	();

	void FindOverlaps(BEDPE &, vector<BED> &, vector<BED> &, vector<BED> &, vector<BED> &, string &);	

	void FindQualityHitsBetweenEnds(BEDPE, int, vector<BED> &, vector<BED> &, int &);
		
private:

	string bedAFilePE;
	string bedBFilePE;
	
	float overlapFraction;
	string searchType;

	// instance of a paired-end bed file class.
	BedFilePE *bedA;

	// instance of a bed file class.
	BedFilePE *bedB;
};

#endif /* PAIRTOPAIR_H */
