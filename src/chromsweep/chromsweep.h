/*****************************************************************************
  chromsweepBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef CHROMSWEEP_H
#define CHROMSWEEP_H

#include "bedFile.h"
// #include "BamReader.h"
// #include "BamWriter.h"
// #include "BamAncillary.h"
// #include "BamAux.h"
// using namespace BamTools;


#include <vector>
#include <queue>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class ChromSweep {

public:

    // constructor
    ChromSweep(string bedAFile, string bedBFile, bool anyHit,
                               bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
                               float overlapFraction, bool noHit, bool writeCount, bool forceStrand,
                               bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput);

    // destructor
    ~ChromSweep(void);
    
    bool Next(pair<BED, vector<BED> > &next);
    
    //pair<BED, vector<BED> > GetNextResult(void);

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

    vector<BED> _cache;
    vector<BED> _hits;
    queue< pair<BED, vector<BED> > > _results;
    
    // variables for the current query and db entries.
    BED _curr_qy, _curr_db;
    BedLineStatus _qy_status, _db_status;
    int _qy_lineNum, _db_lineNum;

    //------------------------------------------------
    // private methods
    //------------------------------------------------
    // void ScanCache(const BED &curr_qy, BedLineStatus qy_status, vector<BED> &db_cache, vector<BED> &hits);
    // 
    // void ChromCheck(BED &curr_qy, BED &curr_db, 
    //                 BedLineStatus &qy_status, BedLineStatus &db_status,
    //                 int &qy_lineNum, int &db_lineNum,
    //                 vector<BED> &db_cache, vector<BED> &hits);

    void ScanCache();
    void ChromCheck();
    
    void Sweep();
    void ReportHits(const BED &curr_qy, const vector<BED> &hits);


};

#endif /* CHROMSWEEP_H */
