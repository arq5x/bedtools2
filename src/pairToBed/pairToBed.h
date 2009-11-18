#ifndef INTERSECTBED_H
#define INTERSECTBED_H

#include "bedFile.h"
#include "bedFilePE.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedIntersectPE {

public:

	// constructor 
	BedIntersectPE(string &, string &, float &, string &);

	// destructor
	~BedIntersectPE(void);

	void FindOverlaps(BEDPE &, vector<BED> &, vector<BED> &, string &); 
	void FindSpanningOverlaps(BEDPE &, vector<BED> &, string &); 
	
	void IntersectBedPE	();
	
	
private:

	string bedAFilePE;
	string bedBFile;
	
	float overlapFraction;
	string searchType;

	// instance of a paired-end bed file class.
	BedFilePE *bedA;

	// instance of a bed file class.
	BedFile *bedB;
};

#endif /* PEINTERSECTBED_H */
