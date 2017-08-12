/* groupBy.h
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#ifndef GROUPBY_H_
#define GROUPBY_H_

#include "ToolBase.h"
#include "ContextGroupBy.h"

class GroupBy : public ToolBase {

public:
    GroupBy(ContextGroupBy *context);
    ~GroupBy();
	virtual bool init(); // after construction
	virtual bool findNext(RecordKeyVector &hits);
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	virtual void cleanupHits(RecordKeyVector &hits);
	//do any last things needed to wrap up.
	virtual bool finalizeCalculations() { return true;}
	virtual void  giveFinalReport(RecordOutputMgr *outputMgr) {}

protected:
	virtual ContextGroupBy *upCast(ContextBase *context) { return static_cast<ContextGroupBy *>(context); }

	vector<int> _groupCols;
	vector<string> _prevFields;
	FileRecordMgr *_queryFRM;
	Record *_prevRecord;
	Record *getNextRecord();
	bool canGroup(Record *);
	void assignPrevFields();
};




#endif /* GROUPBY_H_ */