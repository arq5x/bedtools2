/*
 * mergeFile.h
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#ifndef MERGEFILE_H_
#define MERGEFILE_H_

#include "ToolBase.h"
#include "ContextMerge.h"


class BlockMgr;
class BinTree;

class MergeFile : public ToolBase {

public:
    MergeFile(ContextMerge *context);
    virtual ~MergeFile();
	virtual bool init();
	virtual bool findNext(RecordKeyVector &hits);
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	virtual void cleanupHits(RecordKeyVector &hits);
	virtual bool finalizeCalculations() { return true; }
	virtual void  giveFinalReport(RecordOutputMgr *outputMgr) {}


protected:
	FileRecordMergeMgr *_frm;

	virtual ContextMerge *upCast(ContextBase *context) { return static_cast<ContextMerge *>(context); }


};




#endif /* MERGEFILE_H_ */
