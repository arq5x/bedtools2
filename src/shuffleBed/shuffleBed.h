#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>

using namespace std;

const int MAX_TRIES = 1000000;

//************************************************
// Class methods and elements
//************************************************
class BedShuffle {

public:

// constructor 
	BedShuffle(string &, string &, string &, bool &, bool &, bool &, int &);

// destructor
	~BedShuffle(void);

	void Shuffle();
	void ShuffleWithExclusions();
	
	void ChooseLocus(BED &);

private:

	string bedFile;
	string genomeFile;
	string excludeFile;
	int seed;
	bool sameChrom;
	bool haveExclude;
	bool haveSeed;


	// The BED file from which to compute coverage.
	BedFile *bed;
	BedFile *exclude;

	map<string, int, less<string> > chromSizes;
	vector<string> chroms;
	int numChroms;
};
