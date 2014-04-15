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
#include "NewChromsweep.h"
#include "BinTree.h"
#include "RecordOutputMgr.h"


FileIntersect::FileIntersect(ContextIntersect *context)
: _context(context),
  _blockMgr(NULL),
  _recordOutputMgr(NULL)
{
	_recordOutputMgr = new RecordOutputMgr();
	_recordOutputMgr->init(_context);
	if (_context->getObeySplits()) {
		_blockMgr = new BlockMgr(_context->getOverlapFraction(), _context->getReciprocal());
		_recordOutputMgr->setSplitInfo(_blockMgr);
	}
}

FileIntersect::~FileIntersect(void) {
	delete _blockMgr;
	_blockMgr = NULL;
	delete _recordOutputMgr;
}

void FileIntersect::processHits(RecordKeyList &hits) {
    _recordOutputMgr->printRecord(hits);
}

bool FileIntersect::intersectFiles()
{
	 if (_context->getSortedInput()) {
		 return processSortedFiles();
	 }
	 return processUnsortedFiles();
}

bool FileIntersect::processSortedFiles()
{
    // use the chromsweep algorithm to detect overlaps on the fly.
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
    		processHits(resultSet);
    	} else {
    		processHits(hitSet);
    	}
    }
    return true;
}

bool FileIntersect::processUnsortedFiles()
{
	BinTree *binTree = new BinTree( _context);
	binTree->loadDB();

	FileRecordMgr *queryFRM = _context->getFile(_context->getQueryFileIdx());


	while (!queryFRM->eof()) {
		Record *queryRecord = queryFRM->getNextRecord();
		if (queryRecord == NULL) {
			continue;
		}
		RecordKeyList hitSet(queryRecord);
		binTree->getHits(queryRecord, hitSet);
    	if (_context->getObeySplits()) {
    		RecordKeyList keySet(hitSet.getKey());
    		RecordKeyList resultSet;
    		_blockMgr->findBlockedOverlaps(keySet, hitSet, resultSet);
    		processHits(resultSet);
    	} else {
    		processHits(hitSet);
    	}
		queryFRM->deleteRecord(queryRecord);
	}

	//clean up.
	delete binTree;
	return true;
}
