/*
 * closestHeader.h
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#ifndef CLOSESTHEADER_H_
#define CLOSESTHEADER_H_

#include "intersectFile.h"
#include "ContextClosest.h"
#include "CloseSweep.h"

class ClosestFile : public IntersectFile {

public:
    ClosestFile(ContextClosest *context);
  	bool findNext(RecordKeyVector &hits);
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);

protected:
	virtual ContextClosest *upCast(ContextBase *context) { return static_cast<ContextClosest *>(context); }
	virtual CloseSweep *upCastSweep() { return static_cast<CloseSweep *>(_sweep); }
	virtual void makeSweep();

};




#endif /* CLOSESTHEADER_H_ */
