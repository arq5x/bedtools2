/*
 * summaryFile.h
 *
 *  Created on: Jan 02, 2019
 *      Author: Aaron Quinlan
 */

#ifndef SUMMARYFILE_H_
#define SUMMARYFILE_H_


#include "ToolBase.h"
#include "ContextSummary.h"

struct Interval {
  CHRPOS start;
  CHRPOS end;
};

class SummaryFile : public ToolBase {

public:
    SummaryFile(ContextSummary *context);
    virtual ~SummaryFile();
    virtual bool init();
    virtual bool findNext(RecordKeyVector &hits);
    virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits) {}
    virtual void cleanupHits(RecordKeyVector &hits) {}
    virtual bool finalizeCalculations() { return true; }
    virtual void giveFinalReport(RecordOutputMgr *);


protected:
    FileRecordMgr *_frm;
    Record *_currRec;
    uint64_t _total_length;
    uint64_t _total_intervals;
    const NewGenomeFile *_genomeFile;
    const vector<string> &_chromList;
    map<string, vector<Interval>, std::less<string> > _chromData;

    FileRecordMgr *_inputFile;
    virtual ContextSummary *upCast(ContextBase *context) { return static_cast<ContextSummary *>(context); }

};



#endif /* SUMMARYFILE_H_ */
