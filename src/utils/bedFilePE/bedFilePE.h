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
	
	vector<string> otherFields;
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

	// Open a BEDPE file for reading (creates an istream pointer)
	void Open(void);
	
	// Close an opened BEDPE file.
	void Close(void);
	
	// Get the next BED entry in an opened BED file.
	BedLineStatus GetNextBedPE (BEDPE &bedpe, int &lineNum);
	
	
	// Methods

	void reportBedPETab(const BEDPE &a);
	void reportBedPENewLine(const BEDPE &a);
	void loadBedPEFileIntoMap();
	void splitBedPEIntoBeds(const BEDPE &a, const int &lineNum, BED &bedEntry1, BED &bedEntry2);
		 
		
	void FindOverlapsPerBin(int bEnd, string chrom, int start, int end, string strand, 
		vector<BED> &hits, bool forceStrand);
		
		
	string bedFile;
	unsigned int bedType;
	
	masterBedMap bedMapEnd1;
	masterBedMap bedMapEnd2;
	
private:
	istream *_bedStream;
	
	// methods
	BedLineStatus parseLine (BEDPE &bedpe, const vector<string> &lineVector, int &lineNum);
	bool parseBedPELine (BEDPE &bed, const vector<string> &lineVector, const int &lineNum);
};

#endif /* BEDFILEPE_H */
