/*****************************************************************************
  Jaccard.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "Jaccard.h"
#include "BlockMgr.h"
#include "NewChromsweep.h"

Jaccard::Jaccard(ContextJaccard *context)
: _context(context),
  _intersectionVal(0),
  _unionVal(0),
  _numIntersections(0)
{
	_blockMgr = new BlockMgr(_context->getOverlapFraction(), _context->getReciprocal());
}

Jaccard::~Jaccard(void) {
	delete _blockMgr;
	_blockMgr = NULL;
}

bool Jaccard::calculate() {

	if (!getIntersectionAndUnion()) {
		return false;
	}

	// header
	cout << "intersection\tunion-intersection\tjaccard\tn_intersections" << endl;

	unsigned long adjustedUnion = _unionVal - _intersectionVal;

	cout << _intersectionVal << "\t" << adjustedUnion << "\t" <<
			(float) _intersectionVal / (float)adjustedUnion << "\t" << _numIntersections << endl;
	return true;
}

bool Jaccard::getIntersectionAndUnion() {
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
			_intersectionVal += getTotalIntersection(&resultSet);
		} else {
			_intersectionVal += getTotalIntersection(&hitSet);
		}
	}

	sweep.closeOut();
	unsigned long queryUnion = sweep.getQueryTotalRecordLength();
	unsigned long dbUnion = sweep.getDatabaseTotalRecordLength();

	_unionVal = queryUnion + dbUnion;
	return true;
}

unsigned long Jaccard::getTotalIntersection(RecordKeyList *recList)
{
	unsigned long intersection = 0;
	const Record *key = recList->getKey();
	int keyStart = key->getStartPos();
	int keyEnd = key->getEndPos();

	int hitIdx = 0;
	for (RecordKeyList::const_iterator_type iter = recList->begin(); iter != recList->end(); iter = recList->next()) {
		const Record *currRec = iter->value();
		int maxStart = max(currRec->getStartPos(), keyStart);
		int minEnd = min(currRec->getEndPos(), keyEnd);
		if (_context->getObeySplits()) {
			intersection += _blockMgr->getOverlapBases(hitIdx);
			hitIdx++;
		} else {
			intersection += (unsigned long)(minEnd - maxStart);
		}
	}
	_numIntersections += (int)recList->size();
	return intersection;
}

