/*****************************************************************************
  mapBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "mapBed.h"


// Constructor
BedMap::BedMap(string bedAFile, string bedBFile, int column, string operation,
               float overlapFraction, bool sameStrand, 
               bool diffStrand, bool reciprocal, 
               bool printHeader) {

    _bedAFile            = bedAFile;
    _bedBFile            = bedBFile;
    _column              = column - 1;
    _operation           = operation;
    _overlapFraction     = overlapFraction;
    _sameStrand          = sameStrand;
    _diffStrand          = diffStrand;
    _reciprocal          = reciprocal;
    _printHeader         = printHeader;
    
    Map();
}

// Destructor
BedMap::~BedMap(void) {
}

void BedMap::Map() {

    // create new BED file objects for A and B
    _bedA = new BedFile(_bedAFile);
    _bedB = new BedFile(_bedBFile);

    // use the chromsweep algorithm to detect overlaps on the fly.
    ChromSweep sweep = ChromSweep(_bedB, _bedA, _sameStrand, _diffStrand, _printHeader);

    pair<BED, vector<BED> > hit_set;
    hit_set.second.reserve(100000);
    while (sweep.Next(hit_set)) {
        ApplyHits(hit_set.first, hit_set.second);
    }
}


// ***************************************
// **************   TO DO   **************
// move all of this logic into bedFile.cpp
// so that hits aren't looped through twice.
// ***************************************

void BedMap::ApplyHits(const BED &a, const vector<BED> &hits) {
    // loop through the hits and report those that meet the user's criteria
    vector<BED>::const_iterator h       = hits.begin();
    vector<BED>::const_iterator hitsEnd = hits.end();
    for (; h != hitsEnd; ++h) {
        cout << h->fields[_column] << " " << h->fields.size() << endl;
    }
}



