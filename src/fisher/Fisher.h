#ifndef FISHER_H
#define FISHER_H

#include "bedFile.h"
#include "bedFilePE.h"

#include "ContextFisher.h"

class BlockMgr;

class Fisher {

public:

    Fisher(ContextFisher *context);
    ~Fisher();

    bool calculate();

private:
    ContextFisher *_context;
    BlockMgr *_blockMgr;
    unsigned long _intersectionVal;
    unsigned long _unionVal;
    bool _haveExclude;
    int _numIntersections;
    unsigned long _queryLen;
    unsigned long _dbLen;
    unsigned long _queryCounts;
    unsigned long _dbCounts;
    unsigned long _overlapCounts;
    bool getFisher();
    BedFile *exclude;

    vector<int> _qsizes;
    unsigned long getTotalIntersection(RecordKeyVector &hits);
};

#endif /* FISHER_H */
