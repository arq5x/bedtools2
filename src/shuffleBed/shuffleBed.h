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
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <algorithm>  // for binary search
using namespace std;

const int MAX_TRIES = 1000000;

//************************************************
// Class methods and elements
//************************************************
class BedShuffle {

public:

    // constructor
    BedShuffle(string &bedFile, string &genomeFile, string &excludeFile, string &includeFile, 
                           bool haveSeed, bool haveExclude, bool haveInclude, bool sameChrom, 
                           float overlapFraction, int seed, bool chooseChrom);

    // destructor
    ~BedShuffle(void);

private:

    string _bedFile;
    string _genomeFile;
    string _excludeFile;
    string _includeFile;
    float  _overlapFraction;
    int _seed;
    bool _sameChrom;
    bool _haveExclude;
    bool _haveInclude;
    bool _haveSeed;
    bool _chooseChrom;


    // The BED file from which to compute coverage.
    BedFile *_bed;
    BedFile *_exclude;
    BedFile *_include;

    GenomeFile *_genome;

    vector<string> _chroms;
    int _numChroms;
    vector<string> _includeChroms;
    int _numIncludeChroms;
    uint32_t _genomeSize;

    // methods
    void Shuffle();
    void ShuffleWithExclusions();
    void ShuffleWithInclusions();

    void ChooseLocus(BED &);
    void ChooseLocusFromInclusionFile(BED &);
};
