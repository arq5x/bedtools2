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
#include "Context.h"
#include "FileRecordMgr.h"
#include "NewChromsweep.h"
#include "BinTree.h"
#include "RecordOutputMgr.h"


FileIntersect::FileIntersect(Context *context)
: _context(context),
  _blockMgr(NULL),
  _recordOutputMgr(NULL)
{
	_blockMgr = new BlockMgr();
	_blockMgr->setContext(context);
	_recordOutputMgr = new RecordOutputMgr();
}

FileIntersect::~FileIntersect(void) {
	if (_blockMgr != NULL) {
		delete _blockMgr;
		_blockMgr = NULL;
	}
	delete _recordOutputMgr;
}

void FileIntersect::processHits(RecordKeyList &hits) {
//	if (hits.getKey()->getType() == FileRecordTypeChecker::BAM_RECORD_TYPE) {
//		RecordKeyList blockList(hits.getKey());
//		bool deleteBlocks = false;
//		_blockMgr->getBlocks(blockList, deleteBlocks);
//		_recordOutputMgr->printRecord(hits, &blockList);
//		if (deleteBlocks) {
//			_blockMgr->deleteBlocks(blockList);
//		}
//		return;
//	}
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
    if (!_recordOutputMgr->init(_context)) {
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
	const QuickString &databaseFilename = _context->getDatabaseFileName();
	BinTree *binTree = new BinTree(_context->getDatabaseFileIdx(), _context);

	FileRecordMgr *queryFRM = new FileRecordMgr(_context->getQueryFileIdx(), _context);
	if (!queryFRM->open()) {
		return false;
	}

	if (!binTree->loadDB()) {
		fprintf(stderr, "Error: Unable to load database file %s.\n", databaseFilename.c_str());
		delete binTree;
		exit(1);
	}


    _context->determineOutputType();
    _recordOutputMgr->init(_context);

	while (!queryFRM->eof()) {
		Record *queryRecord = queryFRM->allocateAndGetNextRecord();
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
	queryFRM->close();

	//clean up.
	delete queryFRM;
	delete binTree;
	return true;
}
