#ifndef SUBTRACTBED_H
#define SUBTRACTBED_H

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedSubtract {

public:

	// constructor 
	BedSubtract(string &, string &, float &, bool &);

	// destructor
	~BedSubtract(void);
	
	void reportARemainder(BED &, int &, int &);	
	void reportA(BED &);

	void FindOverlaps(BED &, vector<BED> &); 
	void SubtractBed();
	
	
private:

	string bedAFile;
	string bedBFile;
	float overlapFraction;
	bool noHit;
	bool forceStrand;
	
	// instance of a bed file class.
	BedFile *bedA, *bedB;

};

#endif /* SUBTRACTBED_H */
