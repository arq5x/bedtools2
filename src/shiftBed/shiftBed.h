/*****************************************************************************
  shiftBed.h

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
class BedShift {

public:

    // constructor
    BedShift(string &bedFile, string &genomeFile, float shiftMinus, float shiftPlus, bool fractional, bool printHeader);

    // destructor
    ~BedShift(void);



private:

    string _bedFile;
    string _genomeFile;

    float  _shiftMinus;
    float  _shiftPlus;
    bool   _fractional;
    bool   _printHeader;

    BedFile *_bed;
    GenomeFile *_genome;

    // methods

    void ShiftBed();

    // method to add requested "slop" to a single BED entry
    void AddShift(BED &bed);
};
