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
#include "BamReader.h"
#include "BamWriter.h"
#include "BamAux.h"
using namespace BamTools;

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedIntersect {

public:

	// constructor 
	BedIntersect(string bedAFile, string bedBFile, bool anyHit, 
							   bool writeA, bool writeB, bool writeOverlap, float overlapFraction, 
							   bool noHit, bool writeCount, bool forceStrand, bool reciprocal,
							   bool bamInput, bool bamOutput);

	// destructor
	~BedIntersect(void);
	
	void reportAIntersect(const BED &, int &, int &);	
	void reportA(const BED &);
	void reportB(const BED &);

	bool FindOverlaps(const BED &a, vector<BED> &hits);
	bool FindOneOrMoreOverlap(const BED &a);
	
	void IntersectBed(istream &bedInput);
	void IntersectBam(string bamFile);
		
	void DetermineBedInput();
	
	
private:

	string bedAFile;
	string bedBFile;
	string notInBFile;
	bool anyHit;
	bool writeA;
	bool writeB;
	bool writeCount;
	bool writeOverlap;
	bool forceStrand;
	bool reciprocal;
	float overlapFraction;
	bool noHit;
	bool bamInput;
	bool bamOutput;
	// instance of a bed file class.
	BedFile *bedA, *bedB;

};

#endif /* INTERSECTBED_H */
