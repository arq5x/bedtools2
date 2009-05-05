#ifndef CLOSESTBED_H
#define CLOSESTBED_H

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedClosest {

public:

	// constructor 
	BedClosest(string &, string &, bool &);

	// destructor
	~BedClosest(void);
		
	void reportA(const BED &);
	void reportB(const BED &);
	void reportNullB();

	void ClosestBed();
	void FindWindowOverlaps(BED &, vector<BED> &);
		
private:

	string bedAFile;
	string bedBFile;
	bool forceStrand;
	
	// instance of a bed file class.
	BedFile *bedA, *bedB;

};
#endif /* CLOSEST_H */
