#ifndef BEDFILE_H
#define BEDFILE_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <limits.h>
#include <cstdio>

using namespace std;

//*************************************************
// Genome binning constants
//*************************************************
const int binOffsetsExtended[] =
	{4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};

	
const int _binFirstShift = 17;	    /* How much to shift to get to finest bin. */
const int _binNextShift  = 3;		/* How much to shift to get to next larger bin. */


//*************************************************
// Common data structures
//*************************************************

struct DEPTH {
	unsigned int starts;
	unsigned int ends;
};

/*
	Structure for regular BED3-BED6 records
*/
struct BED {

	// UCSC BED fields
	string chrom;
	int start;
	int end; 
	string name;
	string score;
	string strand;
	
	vector<string> otherFields;

	// Additional fields
	unsigned int count;			// count of number of intervals
								// that overlap this feature
								
	map<unsigned int, DEPTH> depthMap;
	int minOverlapStart;
};



//*************************************************
// Common functions
//*************************************************

int getBin(int, int);
	
// BED sorting
bool sortByChrom(BED const &, BED const &);
bool sortByStart(const BED &, const BED &);
bool byChromThenStart(BED const &, BED const &);

// BED comparsions
int overlaps(const int, const int, const int, const int);
bool leftOf(const int, const int);
int min(const int, int);
int max(const int, int);


// templated function to convert objects to strings
template <typename T>
std::string ToString(const T & value)
{
	std::stringstream ss;
	ss << value;
	return ss.str();
}


// BED Sorting Methods 
bool sortByChrom(BED const &, BED const &);	
bool sortByStart(const BED &, const BED &);
bool sortBySizeAsc(const BED &, const BED &);
bool sortBySizeDesc(const BED &, const BED &);
bool sortByScoreAsc(const BED &, const BED &);
bool sortByScoreDesc(const BED &, const BED &);
bool byChromThenStart(BED const &, BED const &);

//*************************************************
// Common typedefs
//*************************************************

typedef map<int, vector<BED>, std::less<int> > binsToBeds;
typedef map<string, binsToBeds, std::less<string> > masterBedMap;
typedef map<string, vector<BED>, std::less<string> > masterBedMapNoBin;


//************************************************
// BedFile Class methods and elements
//************************************************
class BedFile {

public:

	// Constructor 
	BedFile(string &);

	// Destructor
	~BedFile(void);

	// parse an input line and determine how it should be handled
	bool parseLine (BED &bed, const vector<string> &lineVector, int &lineNum);
	// parse a BED line
	bool parseBedLine (BED &bed, const vector<string> &lineVector, int lineNum);
	// parse a GFF line
	bool parseGffLine (BED &bed, const vector<string> &lineVector, int lineNum);
	
	void loadBedFileIntoMap();
	void loadBedFileIntoMapNoBin();	

	//void binKeeperFind(map<int, vector<BED>, std::less<int> > &, const int, const int, vector<BED> &);
	void binKeeperFind(string chrom, const int start, const int end, vector<BED> &);


	void countHits(map<int, vector<BED>, std::less<int> > &, BED &, bool &);

	// printing methods
	void reportBedTab(const BED &);
	void reportBedNewLine(const BED &);
		
	void reportBedRangeTab(const BED &bed, int start, int end);
	void reportBedRangeNewLine(const BED &bed, int start, int end);
	

	map<string, int> minPosMap;
	map<string, int> maxPosMap;

	// the bedfile with which this instance is associated
	string bedFile;
	unsigned int bedType;  // 3 -6 for BED
						   // 9 for GFF

	vector<BED> bedVector;
	masterBedMap bedMap;
	masterBedMapNoBin bedMapNoBin;
	
private:


	
};

#endif /* BEDFILE_H */
