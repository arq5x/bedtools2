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

	// Additional fields
	unsigned int count;			// count of number of intervals
								// that overlap this feature
								
	map<unsigned int, DEPTH> depthMap;
	int minOverlapStart;
};



//*************************************************
// Common functions
//*************************************************

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
typedef std::map<string, map<int, vector<BED>, std::less<int> >, std::less<string> > masterBedMap;
typedef std::map<string, vector<BED>, std::less<string> > masterBedMapNoBin;


//************************************************
// BedFile Class methods and elements
//************************************************
class BedFile {

public:

	// Constructor 
	BedFile(string &);

	// Destructor
	~BedFile(void);

	// Methods
	bool parseBedLine (BED &, const vector<string> &, int);

	void loadBedFileIntoMap();
	void loadBedFileIntoMapNoBin();	

	void binKeeperFind(map<int, vector<BED>, 
		std::less<int> > &, const int, 
		const int, vector<BED> &);

	void countHits(map<int, vector<BED>, std::less<int> > &, BED &, bool &);

	// printing methods
	void reportBedTab(BED);
	void reportBedNewLine(BED);
		
	void reportBedRangeTab(BED, int, int);
	void reportBedRangeNewLine(BED, int, int);
	
	// a vector of the BED entries in the BED file.
	vector<BED> bedVector;
	masterBedMap bedMap;
	masterBedMapNoBin bedMapNoBin;

	map<string, int> minPosMap;
	map<string, int> maxPosMap;

	// the bedfile with which this instance is associated
	string bedFile;
	unsigned int bedType;

private:
	// none
};

#endif /* BEDFILE_H */
