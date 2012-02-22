/*****************************************************************************
  intersectBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef INTERSECTBED_H
#define INTERSECTBED_H

#include "bedFile.h"
#include "chromsweep.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "BlockedIntervals.h"
#include "BamAncillary.h"
using namespace BamTools;


#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class BedIntersect {

public:

    // constructor
    BedIntersect(string bedAFile, string bedBFile, bool anyHit,
                               bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
                               float overlapFraction, bool noHit, bool leftJoin, bool writeCount, bool sameStrand, bool diffStrand,
                               bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput, bool isUncompressedBam,
                               bool sortedInput, bool printHeader);

    // destructor
    ~BedIntersect(void);

private:

    //------------------------------------------------
    // private attributes
    //------------------------------------------------
    string _bedAFile;
    string _bedBFile;

    bool  _writeA;            // should the original A feature be reported?
    bool  _writeB;            // should the original B feature be reported?
    bool  _writeOverlap;
    bool  _writeAllOverlap;

    bool  _sameStrand;
    bool  _diffStrand;
    bool  _reciprocal;
    float _overlapFraction;

    bool  _anyHit;
    bool  _noHit;
    bool  _leftJoin;
    bool  _writeCount;        // do we want a count of the number of overlaps in B?
    bool  _obeySplits;
    bool  _bamInput;
    bool  _bamOutput;
    bool  _isUncompressedBam;
    bool  _sortedInput;
    bool  _printable;
    bool  _printHeader;
    
    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    //------------------------------------------------
    // private methods
    //------------------------------------------------
    void IntersectBed(istream &bedInput);

    void IntersectBed();

    void IntersectBam(string bamFile);

    bool processHits(const BED &a, const vector<BED> &hits);

    bool FindOverlaps(const BED &a, vector<BED> &hits);

    bool FindBlockedOverlaps(const BED &a, const vector<BED> &a_blocks, 
        const vector<BED> &hits, bool a_is_bam);

    void ReportOverlapDetail(int overlapBases, const BED &a, const BED &b, CHRPOS s, CHRPOS e);

    void ReportOverlapSummary(const BED &a, const int &numOverlapsFound);

};

#endif /* INTERSECTBED_H */
