/*****************************************************************************
  intersectFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "intersectFile.h"
#include "ContextIntersect.h"
#include "FileRecordMgr.h"
#include "BinTree.h"
#include "RecordOutputMgr.h"


IntersectFile::IntersectFile(ContextIntersect *context)
: ToolBase(upCast(context)),
  _sweep(NULL),
  _binTree(NULL),
  _queryFRM(NULL)
{

}

IntersectFile::~IntersectFile(void) {
	delete _sweep;
	_sweep = NULL;

	delete _binTree;
	_binTree = NULL;
}

bool IntersectFile::init() {

	_queryFRM = upCast(_context)->getFile(upCast(_context)->getQueryFileIdx());

	 if (upCast(_context)->getSortedInput()) {
		 makeSweep();
		return _sweep->init();
	 } else {
		_binTree = new BinTree( upCast(_context));
		_binTree->loadDB();
	 }

	 return true;
}

bool IntersectFile::findNext(RecordKeyVector &hits)
{
	 bool retVal = false;
	 if (upCast(_context)->getSortedInput()) {
		retVal = nextSortedFind(hits);
	 } 
	 else {
		retVal = nextUnsortedFind(hits);
	 }

	 if (retVal) {
		 checkSplits(hits);
	 }
	 return retVal;
}

void IntersectFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	outputMgr->printRecord(hits);
}

void IntersectFile::cleanupHits(RecordKeyVector &hits)
{
	_queryFRM->deleteRecord(hits.getKey());
	hits.clearAll();
}

bool IntersectFile::finalizeCalculations()
{
    if (upCast(_context)->getSortedInput() && !upCast(_context)->hasGenomeFile()) 
    {
        if (_context->getNameCheckDisabled())
            _sweep->closeOut(false);
        else
            _sweep->closeOut(true);
    }
    return true;
}

bool IntersectFile::nextSortedFind(RecordKeyVector &hits)
{
    if (!_sweep->next(hits)) {
    	return false;
    }
    return true;
}

bool IntersectFile::nextUnsortedFind(RecordKeyVector &hits)
{

	while (!_queryFRM->eof()) {
		Record *queryRecord = _queryFRM->getNextRecord();
		if (queryRecord == NULL) {
			continue;
		} else {
			_context->testNameConventions(queryRecord);
			hits.setKey(queryRecord);
			 _binTree->getHits(queryRecord, hits);
			return true;
		}
	}
	return false;
}

void IntersectFile::makeSweep() {
	_sweep = new NewChromSweep(upCast(_context));
}

void IntersectFile::checkSplits(RecordKeyVector &hitSet)
{
	if (upCast(_context)->getObeySplits()) {
		RecordKeyVector keySet(hitSet.getKey());
		RecordKeyVector resultSet(hitSet.getKey());
		upCast(_context)->getSplitBlockInfo()->findBlockedOverlaps(keySet, hitSet, resultSet);
		hitSet.swap(resultSet);
	}
}
