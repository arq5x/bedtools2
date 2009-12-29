/*****************************************************************************
  intersectBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
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
	BedIntersect(string &, string &, bool &, bool &, bool &, float &, bool &, bool &, bool &, bool &);

	// destructor
	~BedIntersect(void);
	
	void reportAIntersect(const BED &, int &, int &);	
	void reportA(const BED &);
	void reportB(const BED &);

	void FindOverlaps(BED &, vector<BED> &); 
	void IntersectBed(istream &bedInput);
	void DetermineBedInput();
	
	
private:

	string bedAFile;
	string bedBFile;
	string notInBFile;
	bool anyHit;
	bool writeA;
	bool writeB;
	bool writeCount;
	bool forceStrand;
	bool reciprocal;
	float overlapFraction;
	bool noHit;

	// instance of a bed file class.
	BedFile *bedA, *bedB;

};

#endif /* INTERSECTBED_H */
