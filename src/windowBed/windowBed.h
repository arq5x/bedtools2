/*****************************************************************************
  windowBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
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

	void WindowIntersectBed(istream &bedInput);
	void FindWindowOverlaps(BED &, vector<BED> &);
	void DetermineBedInput();
		
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
