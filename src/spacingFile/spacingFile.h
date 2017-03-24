/*
 * spacingFile.h
 *
 *  Created on: Apr 30, 2015
 *      Author: nek3d
 */

#ifndef SPACINGFILE_H_
#define SPACINGFILE_H_


#include "ToolBase.h"
#include "ContextSpacing.h"

class SpacingFile : public ToolBase {

public:
    SpacingFile(ContextSpacing *context);
    virtual ~SpacingFile();
    virtual bool init();
    virtual bool findNext(RecordKeyVector &hits);
    virtual void processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits);
    virtual void cleanupHits(RecordKeyVector &hits);
    virtual bool finalizeCalculations() { return true; }
    virtual void  giveFinalReport(RecordOutputMgr *) {}


protected:
    Record *_prevRec;
    Record *_currRec;
    string _distance;

    FileRecordMgr *_inputFile;
    virtual ContextSpacing *upCast(ContextBase *context) { return static_cast<ContextSpacing *>(context); }




};



#endif /* SPACINGFILE_H_ */
