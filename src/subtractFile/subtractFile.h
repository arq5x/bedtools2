/*
 * subtractFile.h
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */

#ifndef SUBTRACTFILE_H_
#define SUBTRACTFILE_H_


#include "RecordKeyVector.h"
#include "BlockMgr.h"

using namespace std;

class ContextSubtract;
class BlockMgr;
class RecordOutputMgr;
class BlockMgr;

class SubtractFile {

public:
    SubtractFile(ContextSubtract *context);
    ~SubtractFile(void);

    bool subtractFiles();

private:
    ContextSubtract *_context;
	Record *_queryRec;
	Record *_databaseRec;
	BlockMgr *_blockMgr;
	RecordOutputMgr *_recordOutputMgr;

	void processHits(RecordKeyVector &hits);
	bool processSortedFiles();
	bool processUnsortedFiles();

	void subtractHits(RecordKeyVector &hits);

	BlockMgr *_tmpBlocksMgr;
};




#endif /* SUBTRACTFILE_H_ */
