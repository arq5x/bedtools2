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
#include "BamReader.h"
#include "BamWriter.h"
#include "BamAncillary.h"
#include "BamAux.h"
using namespace BamTools;


#include <vector>
#include <queue>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class BedIntersect {

public:

    // constructor
    BedIntersect(string bedAFile, string bedBFile, bool anyHit,
                               bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
                               float overlapFraction, bool noHit, bool writeCount, bool forceStrand,
                               bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput);

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

    bool  _forceStrand;
    bool  _reciprocal;
    float _overlapFraction;

    bool  _anyHit;
    bool  _noHit;
    bool  _writeCount;        // do we want a count of the number of overlaps in B?
    bool  _obeySplits;
    bool  _bamInput;
    bool  _bamOutput;

    bool _printable;

    queue<BED*> _outputBuffer;
    bool  _lastPick;

    map<string, vector<BED*> > _windowA;
    map<string, vector<BED*> > _windowB;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    //------------------------------------------------
    // private methods
    //------------------------------------------------
    void IntersectBed(istream &bedInput);

    void Scan(BED *x, vector<BED *> *windowX, BedLineStatus xStatus,
        const BED &y, vector<BED *> *windowY, BedLineStatus yStatus);

    void AddHits(BED *x, const BED &y);

    void FlushOutputBuffer(bool final = false);

    vector<BED*>* GetWindow(const string &chrom, bool isA);

    void ChromSwitch(const string &chrom);

    void IntersectBed();

    void IntersectBam(string bamFile);

    bool processHits(BED &a, vector<BED> &hits);

    bool FindOverlaps(const BED &a, vector<BED> &hits);

    bool FindOneOrMoreOverlap(const BED &a);

    void ReportOverlapDetail(const int &overlapBases, const BED &a, const BED &b,
                             const CHRPOS &s, const CHRPOS &e);
    void ReportOverlapSummary(const BED &a, const int &numOverlapsFound);

    void ReportHits(set<BED> &A, set<BED> &B);

};

#endif /* INTERSECTBED_H */
