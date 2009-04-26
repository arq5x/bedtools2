#ifndef CLOSESTBED_H
#define CLOSESTBED_H

#include "bedFile.h"
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <fstream>

using namespace std;


//***********************************************
// Typedefs
//***********************************************
typedef list<BED> bedList;


//************************************************
// Class methods and elements
//************************************************
class BedClosest {

public:

	// constructor 
	BedClosest(string &, string &);

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

	// instance of a bed file class.
	BedFile *bedA, *bedB;

};
#endif /* CLOSEST_H */
