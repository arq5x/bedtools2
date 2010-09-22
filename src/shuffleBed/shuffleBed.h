/*****************************************************************************
  shuffleBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"
#include "genomeFile.h"

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
	BedShuffle(string &bedFile, string &genomeFile, string &excludeFile, 
		bool &haveSeed, bool &haveExclude, bool &sameChrom, int &seed);

	// destructor
	~BedShuffle(void);

private:

	string _bedFile;
	string _genomeFile;
	string _excludeFile;
	int _seed;
	bool _sameChrom;
	bool _haveExclude;
	bool _haveSeed;


	// The BED file from which to compute coverage.
	BedFile *_bed;
	BedFile *_exclude;

	GenomeFile *_genome;

	vector<string> _chroms;
	int _numChroms;
	
	// methods
	void Shuffle();
	void ShuffleWithExclusions();
	
	void ChooseLocus(BED &);
};
