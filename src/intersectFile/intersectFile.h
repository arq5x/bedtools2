/*****************************************************************************
  intersectFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef INTERSECTFILE_H
#define INTERSECTFILE_H

#include "ToolBase.h"
#include "ContextIntersect.h"
#include "NewChromsweep.h"


class BlockMgr;
class BinTree;

class IntersectFile : public ToolBase {

public:
    IntersectFile(ContextIntersect *context);
    virtual ~IntersectFile();
	virtual bool init();
	virtual bool findNext(RecordKeyVector &hits);
	virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
	virtual void cleanupHits(RecordKeyVector &hits);
	virtual bool finalizeCalculations();
	virtual void  giveFinalReport(RecordOutputMgr *outputMgr) {}


protected:
	NewChromSweep *_sweep;
	BinTree *_binTree;
	FileRecordMgr *_queryFRM;

	virtual bool nextSortedFind(RecordKeyVector &hits);
	virtual bool nextUnsortedFind(RecordKeyVector &hits);
	void checkSplits(RecordKeyVector &hits);
	virtual void makeSweep();
	virtual ContextIntersect *upCast(ContextBase *context) { return static_cast<ContextIntersect *>(context); }


};

#endif /* INTERSECTFILE_H */
