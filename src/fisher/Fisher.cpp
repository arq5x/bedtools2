#include "Fisher.h"
#include "BlockMgr.h"
#include "NewChromsweep.h"
#include "kfunc.c"

Fisher::Fisher(ContextFisher *context)
: _context(context),
  _intersectionVal(0),
  _unionVal(0),
  _numIntersections(0),
  _queryLen(0),
  _dbLen(0)
{
    _blockMgr = new BlockMgr(_context->getOverlapFraction(), _context->getReciprocal());
    _haveExclude = false;

    if(!(_context->getExcludeFile().empty())){
        string ex = _context->getExcludeFile();
        exclude = new BedFile(ex);
        exclude->loadBedFileIntoMergedMap();
        _haveExclude = true;
    }
}

Fisher::~Fisher(void) {
    delete _blockMgr;
    _blockMgr = NULL;
}

bool Fisher::calculate() {

    if (!getFisher()) {
        return false;
    }

    // header
    cout << "# Contingency Table" << endl;

    double left, right, two;

    long long genomeSize = _context->getGenomeFile()->getGenomeSize();
    if(_haveExclude){
        genomeSize -= exclude->getTotalFlattenedLength();
    }
    // bases covered by both
    long long n22 = genomeSize - _queryLen - _dbLen + _intersectionVal;
    // bases covered only by -a
    long long n21 = _dbLen - _intersectionVal;
    // bases covered only by -b
    long long n12 = _queryLen - _intersectionVal;
    // bases covered by neither a nor b
    long long n11 = _intersectionVal;

    printf("#_________________________________________\n");
    printf("#           | %-12s | %-12s |\n", " in -b", "not in -b");
    printf("#     in -a | %-12lld | %-12lld |\n", n11, n12);
    printf("# not in -a | %-12lld | %-12lld |\n", n21, n22);
    printf("#_________________________________________\n");

    kt_fisher_exact(n11, n12, n21, n22, &left, &right, &two);
    double ratio = ((double)n11 / (double)n12) / ((double)n21 / (double)n22);

    printf("# p-values for fisher's exact test\n");
    printf("left\tright\ttwo-tail\tratio\n");
    printf("%.5g\t%.5g\t%.5g\t%.3f\n", left, right, two, ratio);
    
    return true;
}

bool Fisher::getFisher() {
    NewChromSweep sweep(_context);
    if (!sweep.init()) {
        return false;
    }

    RecordKeyVector hitSet;
    while (sweep.next(hitSet)) {
        if (_context->getObeySplits()) {
            RecordKeyVector keySet(hitSet.getKey());
            RecordKeyVector resultSet(hitSet.getKey());
            _blockMgr->findBlockedOverlaps(keySet, hitSet, resultSet);
            _intersectionVal += getTotalIntersection(resultSet);
        } else {
            _intersectionVal += getTotalIntersection(hitSet);
        }
    }

    sweep.closeOut();
    _queryLen = sweep.getQueryTotalRecordLength();
    _dbLen = sweep.getDatabaseTotalRecordLength();

    _unionVal = _queryLen + _dbLen;
    return true;
}

unsigned long Fisher::getTotalIntersection(RecordKeyVector &recList)
{
    unsigned long intersection = 0;
    const Record *key = recList.getKey();
    int keyStart = key->getStartPos();
    int keyEnd = key->getEndPos();

    int hitIdx = 0;
    for (RecordKeyVector::const_iterator_type iter = recList.begin(); iter != recList.end(); iter = recList.next()) {
        int maxStart = max((*iter)->getStartPos(), keyStart);
        int minEnd = min((*iter)->getEndPos(), keyEnd);
        if (_context->getObeySplits()) {
            intersection += _blockMgr->getOverlapBases(hitIdx);
            hitIdx++;
        } else {
            intersection += (unsigned long)(minEnd - maxStart);
        }
    }
    _numIntersections += (int)recList.size();
    return intersection;
}

