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
			_recordOutputMgr->printRecord(resultSet.getKey(), _context->getColumnOpsVal(resultSet));
    	} else {
			_recordOutputMgr->printRecord(hitSet.getKey(), _context->getColumnOpsVal(hitSet));
		}
    }
    return true;
}

