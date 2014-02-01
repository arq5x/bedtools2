/*****************************************************************************
  mapFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef MAPFILE_H
#define MAPFILE_H

using namespace std;

#include <sstream>
#include <iomanip>
#include "VectorOps.h"
#include "RecordKeyList.h"

using namespace std;

class ContextMap;
class BlockMgr;
class RecordOutputMgr;

class FileMap {

public:
    FileMap(ContextMap *context);
    ~FileMap(void);

    bool mapFiles();

private:
    ContextMap *_context;
    Record *_queryRec;
    Record *_databaseRec;
    BlockMgr *_blockMgr;
    RecordOutputMgr *_recordOutputMgr;

    vector<string> _column_vec; // vector to hold current column's worth of data

    ostringstream _tmp_output;
    QuickString _output;  // placeholder for the results of mapping B to each a in A.
    //------------------------------------------------
    // private methods
    //------------------------------------------------
    void Map();
    void SummarizeHits(RecordKeyList &hits);
    void ExtractColumnFromHits(RecordKeyList &hits);

};

#endif /* MAPFILE_H */


/*
#include "bedFile.h"
#include "chromsweep.h"
#include "VectorOps.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "BamAncillary.h"
using namespace BamTools;


#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
using namespace std;



class BedMap {

public:

    // constructor
    BedMap(string bedAFile, string bedBFile, int column, string operation,
                   float overlapFraction, bool sameStrand, 
                   bool diffStrand, bool reciprocal, 
                   bool choseNullValue, string nullValue, 
                   bool printHeader);

    // destructor
    ~BedMap(void);

private:

    //------------------------------------------------
    // private attributes
    //------------------------------------------------
    string _bedAFile;
    string _bedBFile;
    int _column;
    string _operation;
    bool  _sameStrand;
    bool  _diffStrand;
    bool  _reciprocal;
    float _overlapFraction;
    string _nullValue;
    bool  _printHeader;
    
    // instance of a bed file class.
    BedFile *_bedA, *_bedB;
    
    vector<string> _column_vec; // vector to hold current column's worth of data

    //------------------------------------------------
    // private methods
    //------------------------------------------------
    void Map();
    string MapHits(const BED &a, const vector<BED> &hits);
    void ExtractColumnFromHits(const vector<BED> &hits);
};
*/
//#endif /* MAPFILE_H */
