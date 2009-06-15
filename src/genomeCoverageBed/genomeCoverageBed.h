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

  void CoverageBeds();

private:
	
	string bedFile;
	string genomeFile;
	bool eachBase;
	bool startSites;
	int max;

	// instance of a bed file class.
	BedFile *bed;
	
	chromDepthMap chromCov;

};
