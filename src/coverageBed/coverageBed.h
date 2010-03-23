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

#include "BamReader.h"
#include "BamAux.h"
using namespace BamTools;

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedCoverage {

public:

	// constructor 
	BedCoverage(string &bedAFile, string &bedBFile, bool &forceStrand, bool &writeHistogram, bool &bamInput);

	// destructor
	~BedCoverage(void);
	
	void CollectCoverageBed(istream &bedInput);

	void CollectCoverageBam(string bamFile);
	
	void DetermineBedInput();
	
private:

	// input files.
	string bedAFile;
	string bedBFile;

	// instance of a bed file class.
	BedFile *bedA, *bedB;
	
	// do we care about strandedness when counting coverage?
	bool forceStrand;
	
	// should we write a histogram for each feature in B?
	bool writeHistogram;
	
	// are we dealing with BAM input for "A"?
	bool bamInput;
	
	// private function for reporting coverage information
	void ReportCoverage();
};
#endif /* COVERAGEBED_H */
