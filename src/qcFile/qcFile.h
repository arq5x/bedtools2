/*
 * qcFile.h
 *
 *  Created on: Jan 02, 2019
 *      Author: Aaron Quinlan
 */

#ifndef QCFILE_H_
#define QCFILE_H_


#include "ToolBase.h"
#include "ContextQc.h"

struct Interval {
  CHRPOS start;
  CHRPOS end;
};

class QcFile : public ToolBase {

public:
    QcFile(ContextQc *context);
    virtual ~QcFile();
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
    virtual ContextQc *upCast(ContextBase *context) { return static_cast<ContextQc *>(context); }

};



#endif /* QCFILE_H_ */
