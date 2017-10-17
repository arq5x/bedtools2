#include "subtractFile.h"
#include <numeric>

SubtractFile::SubtractFile(ContextSubtract *context)
: IntersectFile(context),
  _tmpBlocksMgr(NULL),
  _deleteTmpBlocks(false),
  _dontReport(false)
{
	_tmpBlocksMgr = new BlockMgr(upCast(_context)->getOverlapFractionA(), upCast(_context)->getReciprocalFraction());
	
	// if using -N, we need to set -f to 1E-9 to collect all
	// overlaps, then use what was set as -f for the minimum
	// coverage of A that must be achieved to remove it.
	if (upCast(_context)->getRemoveSum()) {
		upCast(_context)->setSubtractFraction(upCast(_context)->getOverlapFractionA());
    	upCast(_context)->setOverlapFractionA(1E-9);
    }
}

SubtractFile::~SubtractFile() {
	delete _tmpBlocksMgr;
	_tmpBlocksMgr = NULL;
}

bool SubtractFile::findNext(RecordKeyVector &hits)
{
	if (IntersectFile::findNext(hits)) {
		subtractHits(hits);
		return true;
	}
	return false;
}

void SubtractFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	if (!_dontReport) {
		IntersectFile::processHits(outputMgr, hits);
	}
	_dontReport = false;
}

void SubtractFile::cleanupHits(RecordKeyVector &hits)
{
	if (_deleteTmpBlocks) {
	    _tmpBlocksMgr->deleteBlocks(hits);
	    _deleteTmpBlocks = false;
	}
	IntersectFile::cleanupHits(hits);
}

void SubtractFile::subtractHits(RecordKeyVector &hits) {
	if (hits.empty()) {
        // no intersection, nothing to subtract.
        // just copy key to hits as if it were a
        // self-intersection. This is just for reporting
        // purposes.
        hits.push_back(hits.getKey());
		return;
	}

	if (upCast(_context)->getRemoveAll() && upCast(_context)->getSubtractFraction() == 0.0) {
		// hits aren't empty, meaning there is intersection,
		// so we want to not report the hit.
		_dontReport = true;
		return;
	}

	//loop through hits. Track which bases in query were covered
	Record *keyRec = hits.getKey();
	int keyStart = keyRec->getStartPos();
	int keyEnd = keyRec->getEndPos();

	//this vector of bools will represent the bases of the query.
	//for each base, true means uncovered, false means covered.
	//they begin as all uncovered.
	vector<bool> keyBases(keyEnd - keyStart, true);

	//now loop through the hits, and cover corresponding query bases
	//by setting them to false.
	bool basesRemoved = false;
	for (RecordKeyVector::iterator_type iter = hits.begin(); iter != hits.end(); iter = hits.next()) {
		Record *hitRec = *iter;
		int hitStart = hitRec->getStartPos();
		int hitEnd = hitRec->getEndPos();

		int startIdx = max(keyStart, hitStart) - keyStart;
		int endIdx = min(keyEnd, hitEnd) - keyStart;

		int keyLen = keyEnd - keyStart;
		int coveredLen = endIdx - startIdx;
		float coveragePct = (float)coveredLen / (float)keyLen;
		//for each base in the hit, set the base in the query to false.
		//this effectively "erases" the covered bits. Only do
		if (upCast(_context)->getRemoveSum() || coveragePct >= upCast(_context)->getSubtractFraction()) {
			std::fill(keyBases.begin() + startIdx, keyBases.begin() + endIdx, false);
			basesRemoved = true;
		}
	}

	if (!basesRemoved) {
		//treat as if there were no intersection
		hits.clearVector();
		hits.push_back(hits.getKey());
		return;
	} else if (upCast(_context)->getRemoveAll()) {
		_dontReport = true;
		return;
	}
	// if the -N option is used ( removeSum), do not report if the percentage of
	// uniquely covered bases exceeds the overlap fraction.
	if (upCast(_context)->getRemoveSum()) {
		//determine how many bases are left uncovered.
		int numBasesUncovered = std::accumulate(keyBases.begin(), keyBases.end(), 0);
		//determine percentage that are covered.
		float pctCovered = 1.0 - (float)numBasesUncovered / (float)(keyEnd - keyStart);
		if (pctCovered > upCast(_context)->getSubtractFraction()) {
			_dontReport = true;
			return;
		} else {
            hits.clearVector();
            hits.push_back(hits.getKey());
        }
		return;
	}

	//now make "blocks" out of the query's remaining stretches of
	//uncovered bases.
	hits.clearVector();
    for (int i = 0; i < (int)keyBases.size(); i++) {
        if (keyBases[i] == true) {
            int blockStart = keyStart + i;
            while (keyBases[i] == true && i < (int)keyBases.size()) {
                i++;
            }
            int blockEnd = min(keyStart + i, keyEnd);
            hits.push_back(_tmpBlocksMgr->allocateAndAssignRecord(keyRec, blockStart, blockEnd));
        }
    }
    _deleteTmpBlocks = true;

}
