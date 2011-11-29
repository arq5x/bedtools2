/*****************************************************************************
  flankBed.h

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


//************************************************
// Class methods and elements
//************************************************
class BedFlank {

public:

    // constructor
    BedFlank(string &bedFile, string &genomeFile, bool forceStrand, 
             float leftSlop, float rightSlop, bool fractional,
             bool printHeader);

    // destructor
    ~BedFlank(void);



private:

    string _bedFile;
    string _genomeFile;

    bool   _forceStrand;
    float  _leftFlank;
    float  _rightFlank;
    bool   _fractional;
    bool   _printHeader;

    BedFile *_bed;
    GenomeFile *_genome;

    // methods

    void FlankBed();

    // method to grab requested flank w.r.t. a single BED entry
    void AddFlank(BED &bed, int leftSlop, int rightSlop);
    
    // method to grab requested flank w.r.t. a single BED entry, 
    // while choosing flanks based on strand
    void AddStrandedFlank(BED &bed, int leftSlop, int rightSlop);
};
