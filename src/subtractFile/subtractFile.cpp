/*
 * subtractFile.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */


/*****************************************************************************
  subtractFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "subtractFile.h"
#include "ContextSubtract.h"
#include "FileRecordMgr.h"
#include "NewChromsweep.h"
#include "BinTree.h"
#include "RecordOutputMgr.h"

#include <numeric> //for std::accumulate

SubtractFile::SubtractFile(ContextSubtract *context)
: _context(context),
  _blockMgr(NULL),
  _recordOutputMgr(NULL),
  _tmpBlocksMgr(NULL)
{
	_recordOutputMgr = new RecordOutputMgr();
	_recordOutputMgr->init(_context);
	if (_context->getObeySplits()) {
		_blockMgr = new BlockMgr(_context->getOverlapFraction(), _context->getReciprocal());
		_recordOutputMgr->setSplitInfo(_blockMgr);
	}
	_tmpBlocksMgr = new BlockMgr(_context->getOverlapFraction(), _context->getReciprocal());
}

SubtractFile::~SubtractFile(void) {
	delete _blockMgr;
	_blockMgr = NULL;
	delete _recordOutputMgr;
	delete _tmpBlocksMgr;
}

void SubtractFile::processHits(RecordKeyVector &hits) {
    _recordOutputMgr->printRecord(hits);
}

bool SubtractFile::subtractFiles()
{
	 if (_context->getSortedInput()) {
		 return processSortedFiles();
	 }
	 return processUnsortedFiles();
}

bool SubtractFile::processSortedFiles()
{
    // use the chromsweep algorithm to detect overlaps on the fly.
    NewChromSweep sweep(_context);

    if (!sweep.init()) {
    	return false;
    }

    RecordKeyVector hitSet;
    while (sweep.next(hitSet)) {
    	if (_context->getObeySplits()) {
    		RecordKeyVector keySet(hitSet.getKey());
    		RecordKeyVector resultSet(hitSet.getKey());
    		_blockMgr->findBlockedOverlaps(keySet, hitSet, resultSet);
    		subtractHits(resultSet);
    	} else {
    		subtractHits(hitSet);
    	}
    }
    if (!_context->hasGenomeFile()) {
    	sweep.closeOut(true);
    }
    return true;
}

bool SubtractFile::processUnsortedFiles()
{
	BinTree *binTree = new BinTree( _context);
	binTree->loadDB();

	FileRecordMgr *queryFRM = _context->getFile(_context->getQueryFileIdx());


	while (!queryFRM->eof()) {
		Record *queryRecord = queryFRM->getNextRecord();
		if (queryRecord == NULL) {
			continue;
		}
		RecordKeyVector hitSet(queryRecord);
		binTree->getHits(queryRecord, hitSet);
    	if (_context->getObeySplits()) {
    		RecordKeyVector keySet(hitSet.getKey());
    		RecordKeyVector resultSet;
    		_blockMgr->findBlockedOverlaps(keySet, hitSet, resultSet);
    		subtractHits(resultSet);
    	} else {
    		subtractHits(hitSet);
    	}
		queryFRM->deleteRecord(queryRecord);
	}

	//clean up.
	delete binTree;
	return true;
}


void SubtractFile::subtractHits(RecordKeyVector &hits) {
	if (hits.empty()) {
		// no intersection, nothing to subtract.
		// just copy key to hits as if it were a
		// self-intersection. This is just for reporting
		// purposes.
		hits.push_back(hits.getKey());
		processHits(hits);
		return;
	}

	if (_context->getRemoveAll() && _context->getSubtractFraction() == 0.0) {
		// hits aren't empty, meaning there is intersection,
		// so we want to not report the hit.
		hits.clearAll();
		return;
	}

	//loop through hits. Track which bases in query were covered
	const Record *keyRec = hits.getKey();
	int keyStart = keyRec->getStartPos();
	int keyEnd = keyRec->getEndPos();

	//this vector of bools will represent the bases of the query.
	//for each base, true means uncovered, false means covered.
	//they begin as all uncovered.
	vector<bool> keyBases(keyEnd - keyStart, true);

	//now loop through the hits, and cover corresponding query bases
	//by setting them to false.
	bool basesRemoved = false;
	for (RecordKeyVector::const_iterator_type iter = hits.begin(); iter != hits.end(); iter = hits.next()) {
		const Record *hitRec = *iter;
		int hitStart = hitRec->getStartPos();
		int hitEnd = hitRec->getEndPos();

		int startIdx = max(keyStart, hitStart) - keyStart;
		int endIdx = min(keyEnd, hitEnd) - keyStart;

		int keyLen = keyEnd - keyStart;
		int coveredLen = endIdx - startIdx;
		float coveragePct = (float)coveredLen / (float)keyLen;
		//for each base in the hit, set the base in the query to false.
		//this effectively "erases" the covered bits. Only do
		if (_context->getRemoveSum() || coveragePct >= _context->getSubtractFraction()) {
			std::fill(keyBases.begin() + startIdx, keyBases.begin() + endIdx, false);
			basesRemoved = true;
		}
	}

	if (!basesRemoved) {
		//treat as if there were no intersection
		hits.clearVector();
		hits.push_back(hits.getKey());
		processHits(hits);
		return;
	} else if (_context->getRemoveAll()) {
		hits.clearAll();
		return;
	}
	// if the -N option is used ( removeSum), do not report if the percentage of
	// uniquely covered bases exceeds the overlap fraction.
	if (_context->getRemoveSum()) {
		//determine how many bases are left uncovered.
		int numBasesUncovered = std::accumulate(keyBases.begin(), keyBases.end(), 0);
		//determine percentage that are covered.
		float pctCovered = 1.0 - (float)numBasesUncovered / (float)(keyEnd - keyStart);
		if (pctCovered > _context->getSubtractFraction()) {
			hits.clearAll();
			return;
		} else {
			hits.clearVector();
			hits.push_back(hits.getKey());
		}
		processHits(hits);
		return;
	}

	//now make "blocks" out of the query's remaining stretches of
	//uncovered bases.
	RecordKeyVector tempHits(keyRec);
    for (int i = 0; i < (int)keyBases.size(); i++) {
        if (keyBases[i] == true) {
            int blockStart = keyStart + i;
            while (keyBases[i] == true && i < (int)keyBases.size()) {
                i++;
            }
            int blockEnd = min(keyStart + i, keyEnd);
            tempHits.push_back(_tmpBlocksMgr->allocateAndAssignRecord(keyRec, blockStart, blockEnd));
        }
    }
    processHits(tempHits);
    _tmpBlocksMgr->deleteBlocks(tempHits);

}
