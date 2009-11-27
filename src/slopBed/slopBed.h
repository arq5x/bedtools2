#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
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

	void DetermineBedInput();

	void ProcessBed(istream &bedInput);
		
	// method to add requested "slop" to each BED entry
	void AddSlop(BED &);
	
private:

	string bedFile;
	string genomeFile;

	bool forceStrand;
	int leftSlop;
	int rightSlop;
	
	BedFile *bed;
	
	map<string, int, less<string> > chromSizes;
};
