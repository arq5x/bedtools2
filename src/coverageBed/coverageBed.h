/*****************************************************************************
  coverageBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#ifndef	COVERAGEBED_H
#define COVERAGEBED_H

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedCoverage {

public:

	// constructor 
	BedCoverage(string &, string &, bool &);

	// destructor
	~BedCoverage(void);
	
	void GetCoverage(istream &bedInput);
	
	void DetermineBedInput();
	
private:

	string bedAFile;
	string bedBFile;

	// instance of a bed file class.
	BedFile *bedA, *bedB;
	
	bool forceStrand;

};
#endif /* COVERAGEBED_H */
