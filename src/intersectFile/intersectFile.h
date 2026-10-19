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
class WorkerPool;
class SingleLineDelimTextFileReader;

class IntersectFile : public ToolBase {

public:
    IntersectFile(ContextIntersect *context);
    virtual ~IntersectFile();
	virtual bool init();
	virtual bool findNext(RecordKeyVector &);
	virtual void processHits(RecordOutputMgr *, RecordKeyVector &);
	virtual void cleanupHits(RecordKeyVector &);
	virtual bool finalizeCalculations();
	virtual void  giveFinalReport(RecordOutputMgr *) {}
	virtual void setOutputMgr(RecordOutputMgr *outputMgr) { _outputMgr = outputMgr; }


protected:
	NewChromSweep *_sweep;
	BinTree *_binTree;
	FileRecordMgr *_queryFRM;
	// whether processHits gives the output manager each key's line, and
	// the reader whose current line it is
	bool _keyLines;
	const SingleLineDelimTextFileReader *_keyLineReader;
	RecordOutputMgr *_outputMgr;

	// -threads N: its threads, and unsorted queries answered on them
	// (intersectFile.cpp)
	WorkerPool *_pool;
	class ParallelQueries;
	ParallelQueries *_parallel;

	virtual bool nextSortedFind(RecordKeyVector &hits);
	virtual bool nextUnsortedFind(RecordKeyVector &hits);
	virtual void checkSplits(RecordKeyVector &hits);
	virtual void makeSweep();
	virtual ContextIntersect *upCast(ContextBase *context) { return static_cast<ContextIntersect *>(context); }


};

#endif /* INTERSECTFILE_H */
