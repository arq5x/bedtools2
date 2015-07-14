/*
 * fisher.h
 *
 *  Created on: Apr 30, 2015
 *      Author: nek3d
 */

#ifndef FISHER_H_
#define FISHER_H_

#include "jaccard.h"
#include "bedFile.h"
#include "ContextFisher.h"
#include "kfunc.h"

class Fisher : public Jaccard {
public:
    Fisher(ContextFisher *context);
    bool init();
    bool finalizeCalculations();
    void giveFinalReport(RecordOutputMgr *outputMgr);

protected:
    bool _haveExclude;
    unsigned long _queryCounts;
    unsigned long _dbCounts;
    unsigned long _overlapCounts;
    vector<int> _qsizes;
    BedFile *_excludeFile;
    unsigned long getTotalIntersection(RecordKeyVector &hits);
    virtual ContextFisher *upCast(ContextBase *context) { return static_cast<ContextFisher *>(context); }

};



#endif /* FISHER_H_ */
