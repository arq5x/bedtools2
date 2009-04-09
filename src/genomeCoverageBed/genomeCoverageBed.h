#include "bedFile.h"
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <fstream>

using namespace std;


struct COV {
	unsigned short starts;
	unsigned short ends;
	unsigned short depth;
};

//***********************************************
// Typedefs
//***********************************************
typedef list<BED> bedList;


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

};
