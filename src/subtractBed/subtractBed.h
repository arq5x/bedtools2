/*****************************************************************************
  subtractBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
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
	void SubtractBed(istream &bedInput);
	void DetermineBedInput();
	
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
