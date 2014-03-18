/*****************************************************************************
  slopBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "bedFile.h"
#include "GenomeFile.h"

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
    BedSlop(string &bedFile, string &genomeFile, bool forceStrand, 
            float leftSlop, float rightSlop, bool fractional, bool printHeader);

    // destructor
    ~BedSlop(void);



private:

    string _bedFile;
    string _genomeFile;

    bool   _forceStrand;
    float  _leftSlop;
    float  _rightSlop;
    bool   _fractional;
    bool   _printHeader;

    BedFile *_bed;
    GenomeFile *_genome;

    // methods

    void SlopBed();

    // method to add requested "slop" to a single BED entry
    void AddSlop(BED &bed);
};
