#ifndef WINDOWBED_H
#define WINDOWBED_H

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
class BedWindow {

public:

	// constructor 
	BedWindow(string &, string &, int &, bool &, bool &, bool &);

	// destructor
	~BedWindow(void);
		
	void reportA(const BED &);
	void reportB(const BED &);

	void WindowIntersectBed();
	void FindWindowOverlaps(BED &, vector<BED> &);
		
private:

	string bedAFile;
	string bedBFile;
	bool anyHit;
	bool writeCount;
	int slop;
	bool noHit;

	// instance of a bed file class.
	BedFile *bedA, *bedB;

};
#endif /* WINDOWBED_H */
