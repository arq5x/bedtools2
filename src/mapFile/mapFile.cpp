/*****************************************************************************
  mapBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "mapFile.h"
#include "ContextMap.h"
#include "FileRecordMgr.h"
#include "NewChromsweep.h"
#include "BinTree.h"
#include "RecordOutputMgr.h"

const int PRECISION = 21;

FileMap::FileMap(ContextMap *context)
: _context(context),
  _blockMgr(NULL),
  _recordOutputMgr(NULL)
{
  // FIX ME - block manager only works for intersect
  //_blockMgr = new BlockMgr();
  //_blockMgr->setContext(context);
  _recordOutputMgr = new RecordOutputMgr();
}

FileMap::~FileMap(void) {
  if (_blockMgr != NULL) {
    delete _blockMgr;
    _blockMgr = NULL;
  }
  delete _recordOutputMgr;
}

bool FileMap::mapFiles()
{
    NewChromSweep sweep(_context);
    if (!sweep.init()) {
      return false;
    }
    if (!_recordOutputMgr->init(_context)) {
      return false;
    }
    RecordKeyList hitSet;
    while (sweep.next(hitSet)) {
      //if (_context->getObeySplits()) {
      //  RecordKeyList keySet(hitSet.getKey());
      //  RecordKeyList resultSet(hitSet.getKey());
      //  _blockMgr->findBlockedOverlaps(keySet, hitSet, resultSet);
      //} else {
      //}
      //_recordOutputMgr->printKeyAndTerminate(hitSet);
      _recordOutputMgr->printRecord(hitSet.getKey(), MapHits(hitSet));
    }
    return true;
}

void FileMap::ExtractColumnFromHits(RecordKeyList &hits) {
  for (RecordKeyList::const_iterator_type iter = hits.begin(); iter != hits.end(); iter = hits.next()) {
      try {
        ostringstream startString;
        startString << iter->value()->getStartPos();
        //_column_vec.push_back(iter->value().fields.at(_column));
        _column_vec.push_back(startString.str());
      }
      catch(std::out_of_range& e) {
        exit(1);
      }
      /*
      catch(std::out_of_range& e) {
        cerr << endl << "*****" << endl 
             << "*****ERROR: requested column ("
             << _column + 1
             << ") exceeds the number of columns in file at line "
             << _bedA->_lineNum << ". Exiting." 
             << endl << "*****" << endl;
        exit(1);
      }
      */
  }
}

string FileMap::MapHits(RecordKeyList &hits) {
    
    ostringstream output;
    if (hits.size() == 0)
        // FIX ME
        return "-1";
        //return _nullValue;

    ExtractColumnFromHits(hits);

    string operation = _context->getColumnOperation();

    VectorOps vo(_column_vec);
    if (operation == "sum") 
        output << setprecision (PRECISION) << vo.GetSum();
    else if (operation == "mean")
        output << setprecision (PRECISION) << vo.GetMean();
    else if (operation == "median")
        output << setprecision (PRECISION) << vo.GetMedian();
    else if (operation == "min")
        output << setprecision (PRECISION) << vo.GetMin();
    else if (operation == "max")
        output << setprecision (PRECISION) << vo.GetMax();
    else if (operation == "mode")
        output << vo.GetMode();
    else if (operation == "antimode")
        output << vo.GetAntiMode();
    else if (operation == "count") 
        output << setprecision (PRECISION) << vo.GetCount();
    else if (operation == "count_distinct")
        output << setprecision (PRECISION) << vo.GetCountDistinct();
    else if (operation == "collapse")
        output << vo.GetCollapse();
    else if (operation == "distinct")
        output << vo.GetDistinct();
    else {
        cerr << "ERROR: " << operation << " is an unrecoginzed operation\n";
        exit(1);
    }
    _column_vec.clear();
    return output.str();
}


/*
const int PRECISION = 21;
double GetUserColumn(const string s);

// Constructor
BedMap::BedMap(string bedAFile, string bedBFile, int column, string operation,
               float overlapFraction, bool sameStrand, 
               bool diffStrand, bool reciprocal,
               bool choseNullValue, string nullValue, 
               bool printHeader) 
{

    _bedAFile            = bedAFile;
    _bedBFile            = bedBFile;
    _column              = column - 1;  // user's request is 1-based
    _operation           = operation;
    _overlapFraction     = overlapFraction;
    _sameStrand          = sameStrand;
    _diffStrand          = diffStrand;
    _reciprocal          = reciprocal;
    _nullValue           = nullValue;
    _printHeader         = printHeader;
    
    if (!choseNullValue && operation == "count")
        _nullValue = "0";
    Map();
}

// Destructor
BedMap::~BedMap(void) 
{}

void BedMap::Map() {

    // create new BED file objects for A and B
    _bedA = new BedFile(_bedAFile);
    _bedB = new BedFile(_bedBFile);

    // use the chromsweep algorithm to detect overlaps on the fly.
    ChromSweep sweep = ChromSweep(_bedA, _bedB, 
                                  _sameStrand, _diffStrand, 
                                  _overlapFraction, _reciprocal,
                                  false, _printHeader);

    pair<BED, vector<BED> > hit_set;
    hit_set.second.reserve(10000);
    while (sweep.Next(hit_set)) {
        string result = MapHits(hit_set.first, hit_set.second);
        _bedA->reportBedTab(hit_set.first);
        printf("%s\n", result.c_str());
    }
}


string BedMap::MapHits(const BED &a, const vector<BED> &hits) {
    
    ostringstream output;
    if (hits.size() == 0)
        return _nullValue;

    ExtractColumnFromHits(hits);
    VectorOps vo(_column_vec);
    if (_operation == "sum") 
        output << setprecision (PRECISION) << vo.GetSum();
    else if (_operation == "mean")
        output << setprecision (PRECISION) << vo.GetMean();
    else if (_operation == "median")
        output << setprecision (PRECISION) << vo.GetMedian();
    else if (_operation == "min")
        output << setprecision (PRECISION) << vo.GetMin();
    else if (_operation == "max")
        output << setprecision (PRECISION) << vo.GetMax();
    else if (_operation == "mode")
        output << vo.GetMode();
    else if (_operation == "antimode")
        output << vo.GetAntiMode();
    else if (_operation == "count") 
        output << setprecision (PRECISION) << vo.GetCount();
    else if (_operation == "count_distinct")
        output << setprecision (PRECISION) << vo.GetCountDistinct();
    else if (_operation == "collapse")
        output << vo.GetCollapse();
    else if (_operation == "distinct")
        output << vo.GetDistinct();
    else {
        cerr << "ERROR: " << _operation << " is an unrecoginzed operation\n";
        exit(1);
    }
    _column_vec.clear();
    return output.str();
}


void BedMap::ExtractColumnFromHits(const vector<BED> &hits) {
    for (size_t i = 0; i < hits.size(); ++i) {
        try {
            _column_vec.push_back(hits[i].fields.at(_column));
        }
        catch(std::out_of_range& e) {
            cerr << endl << "*****" << endl 
                 << "*****ERROR: requested column ("
                 << _column + 1
                 << ") exceeds the number of columns in file at line "
                 << _bedA->_lineNum << ". Exiting." 
                 << endl << "*****" << endl;
            exit(1);
        }
    }
}
*/
