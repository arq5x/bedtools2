#ifndef WINDOWBED_H
#define WINDOWBED_H

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedWindow {

public:

	// constructor 
	BedWindow(string &, string &, int &, int &, bool &, bool &, bool &, bool &, bool &);

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
	int leftSlop;
	int rightSlop;
	bool noHit;
	bool strandWindows;
	bool matchOnStrand;

	// instance of a bed file class.
	BedFile *bedA, *bedB;

};
#endif /* WINDOWBED_H */
