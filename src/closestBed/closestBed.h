/*****************************************************************************
  closestBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
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
	BedClosest(string &, string &, bool &, string &);

	// destructor
	~BedClosest(void);
		
	void reportA(const BED &);
	void reportB(const BED &);
	void reportNullB();

	void ClosestBed(istream &bedInput);
	void FindWindowOverlaps(BED &, vector<BED> &);
	void DetermineBedInput();
		
private:

	string bedAFile;
	string bedBFile;
	string tieMode;
	bool forceStrand;
	
	// instance of a bed file class.
	BedFile *bedA, *bedB;

};
#endif /* CLOSEST_H */
