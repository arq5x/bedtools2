#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
/*****************************************************************************
  slopBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include <ctime>
using namespace std;


//************************************************
// Class methods and elements
//************************************************
class BedSlop {

public:

	// constructor 
	BedSlop(string &, string &, bool &, int &, int &);

	// destructor
	~BedSlop(void);

	void SlopBed(istream &bedInput);
		
	// method to add requested "slop" to a single BED entry
	void AddSlop(BED &bed);
	
	void DetermineBedInput();
	
private:

	string bedFile;
	string genomeFile;

	bool forceStrand;
	int leftSlop;
	int rightSlop;
	
	BedFile *bed;
	
	map<string, int, less<string> > chromSizes;
};
