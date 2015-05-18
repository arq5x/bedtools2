/*
 * jaccard.h
 *
 *  Created on: Apr 24, 2015
 *      Author: nek3d
 */

#ifndef JACCARD_H_
#define JACCARD_H_

#include "ContextJaccard.h"
#include "intersectFile.h"

class Jaccard : public IntersectFile {

public:
	Jaccard(ContextJaccard *context);
	virtual bool findNext(RecordKeyVector &hits);
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits) {}
	virtual void cleanupHits(RecordKeyVector &hits);
	virtual bool finalizeCalculations();
	virtual void  giveFinalReport(RecordOutputMgr *outputMgr);


protected:
	unsigned long _queryUnion;
	unsigned long _dbUnion;

    unsigned long _intersectionVal;
    unsigned long _unionVal;
    int _numIntersections;

     virtual unsigned long getTotalIntersection(RecordKeyVector &hits);

	virtual ContextJaccard *upCast(ContextBase *context) { return static_cast<ContextJaccard *>(context); }
	virtual FileRecordMergeMgr *upCastFRM(FileRecordMgr *frm) { return static_cast<FileRecordMergeMgr *>(frm); }


};

#endif /* JACCARD_H_ */
