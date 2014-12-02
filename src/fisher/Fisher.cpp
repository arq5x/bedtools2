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
  _dbLen(0),
  _queryCounts(0),
  _dbCounts(0),
  _overlapCounts(0)

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

    double left, right, two;

    long long genomeSize = _context->getGenomeFile()->getGenomeSize();
    if(_haveExclude){
        genomeSize -= exclude->getTotalFlattenedLength();
    }
    // bases covered by neither
    long long n22_full_bases = genomeSize;
    //long long n22_bases = genomeSize - _queryLen - _dbLen + _intersectionVal;
    long double dMean = 1.0 + _dbLen / (long double)_dbCounts;
    long double qMean = 1.0 + _queryLen / (long double)_queryCounts;

    // heursitic, but seems to work quite well -- better than doing more intuitive sum then divide.
    long double bMean = (qMean + dMean);
    //bMean = (_unionVal + 2.0 * _intersectionVal) / (long double)(_dbCounts + _queryCounts);

    long long n11 = (long)_overlapCounts;
    // this could be < 0 because multiple overlaps
    long long n12 = (long)max(0L, (long)_queryCounts - (long)_overlapCounts);
    long long n21 = max(0L, (long)(_dbCounts - _overlapCounts));
    long long n22_full = max(n21 + n21 + n11, (long long)(n22_full_bases / bMean));
    long long n22 = max(0L, (long)(n22_full - n12 - n21 - n11));

    printf("# Number of query intervals: %lu\n", _queryCounts);
    printf("# Number of db intervals: %lu\n", _dbCounts);
    printf("# Number of overlaps: %lu\n", _overlapCounts);
    printf("# Number of possible intervals (estimated): %lld\n", n22_full);

    printf("# phyper(%lld - 1, %lu, %lld - %lu, %lu, lower.tail=F)\n", n11, _queryCounts, n22_full, _queryCounts, _dbCounts);
    cout << "# Contingency Table Of Counts" << endl;
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

    _queryCounts = sweep.getQueryTotalRecords();
    _dbCounts = sweep.getDatabaseTotalRecords();

    _unionVal = _queryLen + _dbLen;
    return true;
}

unsigned long Fisher::getTotalIntersection(RecordKeyVector &recList)
{
    unsigned long intersection = 0;
    const Record *key = recList.getKey();
    int keyStart = key->getStartPos();
    int keyEnd = key->getEndPos();

    _overlapCounts += recList.size();
    _qsizes.push_back((keyEnd - keyStart));

    int hitIdx = 0;
    for (RecordKeyVector::const_iterator_type iter = recList.begin(); iter != recList.end(); iter = recList.next()) {
        int maxStart = max((*iter)->getStartPos(), keyStart);
        int minEnd = min((*iter)->getEndPos(), keyEnd);
        _qsizes.push_back((int)(minEnd - maxStart));
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
