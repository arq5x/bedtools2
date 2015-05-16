/*
 * mergeFile.cpp
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#include "mergeFile.h"

MergeFile::MergeFile(ContextMerge *context)
: ToolBase(context) {

}

MergeFile::~MergeFile() {

}

bool MergeFile::init()
{
	_frm = static_cast<FileRecordMergeMgr *>(upCast(_context)->getFile(0));
	return true;
}

bool MergeFile::findNext(RecordKeyVector &hits)
{
    while (!_frm->eof()) {
    	_frm->getNextRecord(&hits);
    	if (hits.getKey() == NULL) continue;
    	return true;
    }
    return false;

}

void MergeFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	outputMgr->printRecord(hits.getKey(), upCast(_context)->getColumnOpsVal(hits));

}

void MergeFile::cleanupHits(RecordKeyVector &hits)
{
	_frm->deleteMergedRecord(hits);
}

