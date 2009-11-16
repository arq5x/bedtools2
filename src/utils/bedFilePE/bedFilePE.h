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
#include "bedFile.h"
#include "lineFileUtilities.h"

using namespace std;


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
	string score;
	
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
	void reportBedPETab(const BEDPE &);
	void reportBedPENewLine(const BEDPE &);
	void loadBedPEFileIntoMap();
	void splitBedPEIntoBeds(const BEDPE &, unsigned int, BED &, BED &);
		 
	void binKeeperFind(map<int, vector<BED>, 
		std::less<int> > &, const int, 
		const int, vector<BED> &);	
		
	string bedFile;
	unsigned int bedType;
	
	masterBedMap bedMapEnd1;
	masterBedMap bedMapEnd2;
	
private:
	// none
};

#endif /* BEDFILEPE_H */
