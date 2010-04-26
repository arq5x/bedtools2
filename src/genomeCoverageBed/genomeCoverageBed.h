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
#include "genomeFile.h"

#include "BamReader.h"
#include "BamAux.h"
using namespace BamTools;

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;


//***********************************************
// Typedefs
//***********************************************
typedef map<int, DEPTH, less<int> >          depthMap;
typedef map<string, depthMap, less<string> > chromDepthMap;

typedef map<int, unsigned int, less<int> >   histMap;
typedef map<string, histMap, less<string> >  chromHistMap;

//************************************************
// Class methods and elements
//************************************************
class BedGenomeCoverage {

public:

	// constructor 
	BedGenomeCoverage(string bedFile, string genomeFile, bool eachBase, bool startSites, 
		bool bedGraph, int max, bool bamInput);

	// destructor
	~BedGenomeCoverage(void);

private:

	// data
	string _bedFile;
	string _genomeFile;
	bool   _bamInput;
	bool   _eachBase;
	bool   _startSites;
	bool   _bedGraph;
	int    _max;

	// The BED file from which to compute coverage.
	BedFile    *_bed;
	GenomeFile *_genome;
	
	chromDepthMap _chromCov;
	
	// methods
	void CoverageBed();
	void CoverageBam(string bamFile);
	void ReportChromCoverage(const vector<DEPTH> &, const int &chromSize, const string &chrom, chromHistMap&);
	void ReportGenomeCoverage(chromHistMap &chromDepthHist);
	void ReportChromCoverageBedGraph(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom);
	
};
