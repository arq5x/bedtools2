#ifndef BEDFILEPE_H
#define BEDFILEPE_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>

using namespace std;

//*************************************************
// Common data structures
//*************************************************

/*
	Structure for paired-end records
*/
struct BEDPE {

	// UCSC BED fields
	string chrom1;
	int start1;
	int end1;
	
	string chrom2;
	int start2;
	int end2;
	 
	string name;
	int score;
	
	string strand1;
	string strand2;
};



//************************************************
// BedFile Class methods and elements
//************************************************
class BedFilePE {

public:

	// Constructor 
	BedFilePE(string &);

	// Destructor
	~BedFilePE(void);

	// Methods
	bool parseBedPELine (BEDPE &, const vector<string> &, const int &);
	void reportBedPE(const BEDPE &);
	 
	string bedFile;
	unsigned int bedType;

private:
	// none
};

#endif /* BEDFILEPE_H */
