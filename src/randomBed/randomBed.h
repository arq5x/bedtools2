/*****************************************************************************
  randomBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
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
class BedRandom {

public:

    // constructor
    BedRandom(string &genomeFile, uint32_t numToGenerate, int seed,
               bool haveSeed, uint32_t length);

    // destructor
    ~BedRandom(void);

private:

    string _genomeFile;
    int _seed;
    bool _haveSeed;

    GenomeFile *_genome;
    uint32_t _length;
    uint32_t _numToGenerate;
    
    // methods
    void Generate();

};
