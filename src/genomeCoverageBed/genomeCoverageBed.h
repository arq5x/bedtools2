/*****************************************************************************
genomeCoverage.h

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;


//***********************************************
// Typedefs
//***********************************************
typedef map<int, DEPTH, less<int> > depthMap;
typedef map<string, depthMap, less<string> > chromDepthMap;

typedef map<int, unsigned int, less<int> > histMap;
typedef map<string, histMap, less<string> > chromHistMap;

//************************************************
// Class methods and elements
//************************************************
class BedCoverage {

public:

	// constructor 
	BedCoverage(string &, string &, bool &, bool &, int &);

	// destructor
	~BedCoverage(void);

	void CoverageBeds(istream &bedInput);
	void ReportChromCoverage(vector<DEPTH> &, int &, string &, chromHistMap&);
	void ReportGenomeCoverage(map<string, int> &, chromHistMap&);

	void DetermineBedInput();

private:

	string bedFile;
	string genomeFile;
	bool eachBase;
	bool startSites;
	int max;

	// The BED file from which to compute coverage.
	BedFile *bed;
	chromDepthMap chromCov;

};
