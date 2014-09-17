/*****************************************************************************
  NewChromsweep.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef NEW_CHROMSWEEP_H
#define NEW_CHROMSWEEP_H

#include <string>
#include "BTlist.h"
#include "RecordKeyList.h"
#include "RecordKeyVector.h"
#include <queue>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "QuickString.h"

class Record;
class FileRecordMgr;
class ContextIntersect;

class NewChromSweep {
public:

    NewChromSweep(ContextIntersect *context, bool useMergedIntervals = false);
    
    
    ~NewChromSweep(void);
    bool init();
    
    typedef RecordList recListType;
    typedef const RecordListNode *recListIterType;
    // loads next (a pair) with the current query and it's overlaps
    bool next(RecordKeyVector &next);
    
    // NOTE! You MUST call this method after sweep if you want the
    // getTotalRecordLength methods to return the whole length of the
    // record files, rather than just what was used by sweep.
    void closeOut();


    unsigned long getQueryTotalRecordLength() { return _queryRecordsTotalLength; }
    unsigned long getDatabaseTotalRecordLength() { return _databaseRecordsTotalLength; }

private:
    ContextIntersect *_context;
    FileRecordMgr *_queryFRM;
    int _numDBs; //don't really need this stored, but here for code brevity.
    vector<FileRecordMgr *> _dbFRMs;

    bool _useMergedIntervals;

    unsigned long _queryRecordsTotalLength;
    vector<unsigned long> _dbFileRecordsLength; //each value in this vector have the
    //length of all records in the corresponding db file.

    unsigned long _databaseRecordsTotalLength;

    bool _wasInitialized;

    // a cache of still active features from the database file
//    typedef enum { BEFORE_QUERY, NEAR_QUERY, AFTER_QUERY, EOF } cacheStatusType;
//    vector <pair<cacheStatusType, recListType> >_caches;

    vector <recListType>_caches;
    // the set of hits in the database for the current query
//    recListType _hits;

    // the current query and db features.
    const Record * _currQueryRec;
    vector<const Record *> _currDbRecs;

    // a cache of the current chrom from the query. used to handle chrom changes.
    QuickString _currChromName;
    bool _runToQueryEnd;

    void masterScan(RecordKeyVector &retList);

    bool nextRecord(bool query, int dbIdx = -1); //true fetches next query record, false fetches next db record.
    
    void scanCache(int dbIdx, RecordKeyVector &retList);
    void clearCache(int dbIdx);
    bool chromChange(int dbIdx, RecordKeyVector &retList);

    bool dbFinished(int dbIdx);

    bool intersects(const Record *rec1, const Record *rec2) const;

    bool allCachesEmpty();
    bool allCurrDBrecsNull();

};

#endif /* NewChromSweep_H */
