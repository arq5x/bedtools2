/*
 * ToolBase.h
 *
 *  Created on: Mar 25, 2015
 *      Author: nek3d
 */

#ifndef TOOLBASE_H_
#define TOOLBASE_H_

#include "RecordKeyVector.h"
#include "RecordOutputMgr.h"

using namespace std;

class ContextBase;

class ToolBase {
public:
	ToolBase(ContextBase *context) {_context = context; }
	virtual ~ToolBase() {}
	virtual bool init() = 0; // after construction
	virtual bool findNext(RecordKeyVector &hits) = 0;
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits) = 0;
	virtual void cleanupHits(RecordKeyVector &hits) = 0;
	//do any last things needed to wrap up.
	virtual bool finalizeCalculations() = 0;
	virtual void  giveFinalReport(RecordOutputMgr *outputMgr) = 0;
protected:
	ContextBase *_context;
};



#endif /* TOOLBASE_H_ */
