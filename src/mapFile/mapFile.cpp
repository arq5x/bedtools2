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
  _blockMgr = new BlockMgr(_context->getOverlapFraction(), _context->getReciprocal());
  _recordOutputMgr = new RecordOutputMgr();
  _recordOutputMgr->init(_context);
}

FileMap::~FileMap(void) {
	delete _blockMgr;
	_blockMgr = NULL;
	delete _recordOutputMgr;
	_recordOutputMgr = NULL;
}

bool FileMap::mapFiles()
{
    NewChromSweep sweep(_context);
    if (!sweep.init()) {
      return false;
    }
    RecordKeyList hitSet;
    while (sweep.next(hitSet)) {
    	if (_context->getObeySplits()) {
			RecordKeyList keySet(hitSet.getKey());
			RecordKeyList resultSet(hitSet.getKey());
			_blockMgr->findBlockedOverlaps(keySet, hitSet, resultSet);
			SummarizeHits(resultSet);
			_recordOutputMgr->printRecord(resultSet.getKey(), _output);
    	} else {
			SummarizeHits(hitSet);
			_recordOutputMgr->printRecord(hitSet.getKey(), _output);
		}
    }
    return true;
}

void FileMap::ExtractColumnFromHits(RecordKeyList &hits) {
  _column_vec.clear();
  RecordKeyList::const_iterator_type iter = hits.begin();
  for (; iter != hits.end(); iter = hits.next()) 
  {
    _column_vec.push_back(iter->value()->getField(_context->getColumn()).str());
  }
} 

void FileMap::SummarizeHits(RecordKeyList &hits) {

    const QuickString & operation = _context->getColumnOperation();
    _output.clear();

    if (hits.size() == 0) {
        if (operation == "count" || operation == "count_distinct")
            _output.append("0");
        else
            _output.append(_context->getNullValue().str());
        return;
    } 

    _tmp_output.str("");
    _tmp_output.clear();

    ExtractColumnFromHits(hits);

    VectorOps vo(_column_vec);
    if (operation == "sum") 
        _tmp_output << setprecision (PRECISION) << vo.GetSum();
    else if (operation == "mean")
        _tmp_output << setprecision (PRECISION) << vo.GetMean();
    else if (operation == "median")
        _tmp_output << setprecision (PRECISION) << vo.GetMedian();
    else if (operation == "min")
        _tmp_output << setprecision (PRECISION) << vo.GetMin();
    else if (operation == "max")
        _tmp_output << setprecision (PRECISION) << vo.GetMax();
    else if (operation == "absmin")
        _tmp_output << setprecision (PRECISION) << vo.GetAbsMin();
    else if (operation == "absmax")
        _tmp_output << setprecision (PRECISION) << vo.GetAbsMax();
    else if (operation == "mode")
        _tmp_output << vo.GetMode();
    else if (operation == "antimode")
        _tmp_output << vo.GetAntiMode();
    else if (operation == "count") 
        _tmp_output << setprecision (PRECISION) << vo.GetCount();
    else if (operation == "count_distinct")
        _tmp_output << setprecision (PRECISION) << vo.GetCountDistinct();
    else if (operation == "collapse")
        _tmp_output << vo.GetCollapse();
    else if (operation == "distinct")
        _tmp_output << vo.GetDistinct();
    else {
        cerr << "ERROR: " << operation << " is an unrecognized operation\n";
        exit(1);
    }
    _output.append(_tmp_output.str());

}
