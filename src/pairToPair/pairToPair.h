/*****************************************************************************
  pairToPair.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef PAIRTOPAIR_H
#define PAIRTOPAIR_H

#include "bedFile.h"
#include "bedFilePE.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;



//************************************************
// Class methods and elements
//************************************************
class PairToPair {

public:

    // constructor
    PairToPair(string &bedAFilePE, string &bedBFilePE, float &overlapFraction,
        string searchType, bool ignoreStrand, bool reqDiffNames, int slop, bool strandedSlop);

    // destructor
    ~PairToPair(void);

    void IntersectPairs();


private:

    string _bedAFilePE;
    string _bedBFilePE;

    float _overlapFraction;
    string _searchType;
    bool _ignoreStrand;
    bool _reqDiffNames;
    int _slop;
    bool _strandedSlop;

    // instance of a paired-end bed file class.
    BedFilePE *_bedA;

    // instance of a bed file class.
    BedFilePE *_bedB;

    // methods
    // void FindOverlaps(const BEDPE &a, vector<MATE> &hitsA1B1, vector<MATE> &hitsA1B2,
    //  vector<MATE> &hitsA2B1, vector<MATE> &hitsA2B2);
    void FindOverlaps(const BEDPE &a);

    void FindQualityHitsBetweenEnds(CHRPOS start, CHRPOS end,
        const vector<MATE> &hits, vector<MATE> &qualityHits, int &numOverlaps);

    bool FindHitsOnBothEnds(const BEDPE &a, const vector<MATE> &qualityHitsEnd1,
        const vector<MATE> &qualityHitsEnd2);

    void FindHitsOnEitherEnd(const BEDPE &a, const vector<MATE> &qualityHitsEnd1,
        const vector<MATE> &qualityHitsEnd2);

};

#endif /* PAIRTOPAIR_H */
