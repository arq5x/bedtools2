/*
 * subtractFile.h
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */

#ifndef SUBTRACTFILE_H_
#define SUBTRACTFILE_H_


#include "intersectFile.h"
#include "ContextSubtract.h"

class SubtractFile : public IntersectFile {

public:
    SubtractFile(ContextSubtract *context);
    virtual ~SubtractFile();

	virtual bool findNext(RecordKeyVector &hits);
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	virtual void cleanupHits(RecordKeyVector &hits);

protected:
	BlockMgr *_tmpBlocksMgr;
	bool _deleteTmpBlocks;
	bool _dontReport;

	virtual ContextSubtract *upCast(ContextBase *context) { return static_cast<ContextSubtract *>(context); }
	void subtractHits(RecordKeyVector &hits);
};



#endif /* SUBTRACTFILE_H_ */
