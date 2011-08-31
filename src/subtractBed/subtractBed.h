/*****************************************************************************
  subtractBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef SUBTRACTBED_H
#define SUBTRACTBED_H

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedSubtract {

public:

    // constructor
    BedSubtract(string &bedAFile, string &bedBFile, float overlapFraction, bool sameStrand, bool diffStrand);

    // destructor
    ~BedSubtract(void);

private:

    // processing variables
    string _bedAFile;
    string _bedBFile;
    float _overlapFraction;
    bool _sameStrand;
    bool _diffStrand;


    // instances of bed file class.
    BedFile *_bedA, *_bedB;

    // methods
    void FindAndSubtractOverlaps(BED &a, vector<BED> &hits);
    void SubtractBed();
};

#endif /* SUBTRACTBED_H */
