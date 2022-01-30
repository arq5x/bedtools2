/*
 * fisher.cpp
 *
 *  Created on: Apr 30, 2015
 *      Author: nek3d
 */

#include "fisher.h"

Fisher::Fisher(ContextFisher *context)
: Jaccard(context),
  _haveExclude(false),
  _queryCounts(0),
  _dbCounts(0),
  _overlapCounts(0),
  _excludeFile(NULL)

{
	context->runToDbEnd();
}

bool Fisher::init(void)
{
	if(!(upCast(_context)->getExcludeFile().empty())){
		string ex = upCast(_context)->getExcludeFile();
		_excludeFile = new BedFile(ex);
		_excludeFile->loadBedFileIntoMergedMap();
		_haveExclude = true;
	}
	return Jaccard::init();
}

bool Fisher::finalizeCalculations()
{
    _sweep->closeOut();
    _queryUnion = _sweep->getQueryTotalRecordLength();
    _dbUnion = _sweep->getDatabaseTotalRecordLength();

    _queryCounts = _sweep->getQueryTotalRecords();
    _dbCounts = _sweep->getDatabaseTotalRecords();

    _unionVal = _queryUnion + _dbUnion;
    return true;
}

void Fisher::giveFinalReport(RecordOutputMgr *outputMgr)
{
    double left, right, two;

    long long genomeSize = _context->getGenomeFile()->getGenomeSize();
    if(_haveExclude){
        genomeSize -= _excludeFile->getTotalFlattenedLength();
    }
    // bases covered by neither
    long long n22_full_bases = genomeSize;
    //long long n22_bases = genomeSize - _queryUnion - _dbUnion + _intersectionVal;
    long double dMean = 1.0 + _dbUnion / (long double)_dbCounts;
    long double qMean = 1.0 + _queryUnion / (long double)_queryCounts;

    // heursitic, but seems to work quite well -- better than doing more intuitive sum then divide.
    long double bMean = (qMean + dMean);
    //bMean = (_unionVal + 2.0 * _intersectionVal) / (long double)(_dbCounts + _queryCounts);

    long long n11 = (long)_overlapCounts;
    // this could be < 0 because multiple overlaps
    long long n12 = (long)max(0L, (long)_queryCounts - (long)_overlapCounts);
    long long n21 = max(0L, (long)(_dbCounts - _overlapCounts));
    long long n22_full = max(n21 + n12 + n11, (long long)(n22_full_bases / bMean));
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
    
    /* Some implementations report NAN as negative, some as positive. To ensure
     * we get consistent output from each compiler, we should do our own test.
     * Since the test script assumes "-nan", let's setle on that.
     */
    if(std::isnan(ratio)) {
        printf("%.5g\t%.5g\t%.5g\t-nan\n", left, right, two);
    } else {
        printf("%.5g\t%.5g\t%.5g\t%.3f\n", left, right, two, ratio);
    }
}

unsigned long Fisher::getTotalIntersection(RecordKeyVector &recList)
{
    unsigned long intersection = 0;
    Record *key = recList.getKey();
    CHRPOS keyStart = key->getStartPos();
    CHRPOS keyEnd = key->getEndPos();

    _overlapCounts += recList.size();
    // note that we truncate to a max size of 2.1GB
    _qsizes.push_back((int)(keyEnd - keyStart));

    int hitIdx = 0;
    for (RecordKeyVector::iterator_type iter = recList.begin(); iter != recList.end(); iter = recList.next()) {
        CHRPOS maxStart = max((*iter)->getStartPos(), keyStart);
        CHRPOS minEnd = min((*iter)->getEndPos(), keyEnd);
        _qsizes.push_back((int)(minEnd - maxStart));
        if (_context->getObeySplits()) {
            intersection += upCast(_context)->getSplitBlockInfo()->getOverlapBases(hitIdx);
            hitIdx++;
        } else {
            intersection += (unsigned long)(minEnd - maxStart);
        }
    }
    _numIntersections += (int)recList.size();
    return intersection;
}
