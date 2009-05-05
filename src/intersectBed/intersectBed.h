#ifndef INTERSECTBED_H
#define INTERSECTBED_H

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedIntersect {

public:

	// constructor 
	BedIntersect(string &, string &, bool &, bool &, bool &, float &, bool &, bool &, bool &);

	// destructor
	~BedIntersect(void);
	
	void reportAIntersect(const BED &, int &, int &);	
	void reportA(const BED &);
	void reportB(const BED &);

	void FindOverlaps(BED &, vector<BED> &); 
	void IntersectBed();
	
	
private:

	string bedAFile;
	string bedBFile;
	string notInBFile;
	bool anyHit;
	bool writeA;
	bool writeB;
	bool writeCount;
	bool forceStrand;
	float overlapFraction;
	bool noHit;

	// instance of a bed file class.
	BedFile *bedA, *bedB;

};

#endif /* INTERSECTBED_H */
