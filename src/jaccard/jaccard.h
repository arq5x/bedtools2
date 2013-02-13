/*****************************************************************************
  jaccard.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef JACCARD_H
#define JACCARD_H

#include "bedFile.h"
#include "chromsweep.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BlockedIntervals.h"
#include "BamAncillary.h"
using namespace BamTools;


#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class Jaccard {

public:

    // constructor
    Jaccard(string bedAFile, string bedBFile, 
            float overlapFraction, bool reciprocal);

    // destructor
    ~Jaccard(void);

private:

    //------------------------------------------------
    // private attributes
    //------------------------------------------------
    string _bedAFile;
    string _bedBFile;

    bool  _reciprocal;
    float _overlapFraction;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    //------------------------------------------------
    // private methods
    //------------------------------------------------
    unsigned long GetIntersection(size_t &n_intersections);
    unsigned long GetUnion();
    void CalculateJaccard();

    size_t GetTotalIntersection(const BED &a, const vector<BED> &hits);

};

#endif /* JACCARD_H */
