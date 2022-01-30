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
#include "string.h"

using namespace std;

class Record;
class FileRecordMgr;
class ContextIntersect;

class NewChromSweep {
public:

    NewChromSweep(ContextIntersect *context);
    
    
    virtual ~NewChromSweep(void);
    virtual bool init();
    
    typedef RecordList recListType;
    typedef const RecordListNode *recListIterType;

    virtual bool next(RecordKeyVector &next);
    
    // NOTE! You MUST call this method after sweep if you want the
    // getTotalRecordLength methods to return the whole length of the
    // record files, rather than just what was used by sweep.
    void closeOut(bool testChromOrder = false);


    unsigned long getQueryTotalRecordLength() { return _queryRecordsTotalLength; }
    unsigned long getDatabaseTotalRecordLength() { return _databaseRecordsTotalLength; }

    unsigned long getQueryTotalRecords() { return _queryTotalRecords; }
    unsigned long getDatabaseTotalRecords() { return _databaseTotalRecords; }


protected:
    ContextIntersect *_context;
    FileRecordMgr *_queryFRM;
    int _numDBs; //don't really need this stored, but here for code brevity.
    int _numFiles; //ditto. Just numDBs + num queries, which for now is always 1.
    vector<FileRecordMgr *> _dbFRMs;

     unsigned long _queryRecordsTotalLength;
    vector<unsigned long> _dbFileRecordsLength; //each value in this vector have the
    //length of all records in the corresponding db file.

    unsigned long _databaseRecordsTotalLength;

    unsigned long _queryTotalRecords;
    unsigned long _databaseTotalRecords;


    bool _wasInitialized;

    // a cache of still active features from the database file
//    typedef enum { BEFORE_QUERY, NEAR_QUERY, AFTER_QUERY, EOF } cacheStatusType;
//    vector <pair<cacheStatusType, recListType> >_caches;

    vector <recListType>_caches;
    // the set of hits in the database for the current query
//    recListType _hits;

    // the current query and db features.
    Record * _currQueryRec;
    vector<Record *> _currDbRecs;

    // a cache of the current chrom from the query. used to handle chrom changes.
    string _currQueryChromName;
    string _prevQueryChromName;
    bool _runToQueryEnd;
	bool _runToDbEnd;


    virtual void masterScan(RecordKeyVector &retList);

    bool nextRecord(bool query, int dbIdx = -1); //true fetches next query record, false fetches next db record.
    
    virtual void scanCache(int dbIdx, RecordKeyVector &retList);
    virtual void clearCache(int dbIdx);
    virtual bool chromChange(int dbIdx, RecordKeyVector &retList, bool wantScan);

    bool dbFinished(int dbIdx);

    bool intersects(const Record *rec1, const Record *rec2) const;

    bool allCachesEmpty();
    bool allCurrDBrecsNull();


    //
    // members and methods for detecting differently
    // sorted files without a genome file.
    //

    typedef map<string, int> _orderTrackType;
    vector<_orderTrackType *> _fileTracks;
    vector<char*> _filePrevChrom;
    bool _lexicoDisproven; //whether we've established that any file ISN'T in lexicographical order
    bool _lexicoAssumed; //whether we've had to try to guess that any file might be in lexicographical order.
    string _lexicoAssumedChromName; //which chromosome we had to make that guess for. Used in error reporting.
    int _lexicoAssumedFileIdx; //which file we had to make the guess for. Also for error reporting.
    bool _testLastQueryRec;

    void testChromOrder(const Record *rec);
    bool queryChromAfterDbRec(const Record *dbRec);
    int findChromOrder(const Record *rec);
    bool verifyChromOrderMismatch(const string & chrom, const string &prevChrom, int skipFile);
    void testThatAllDbChromsExistInQuery();
    bool testLexicoQueryAfterDb(const Record *queryRec, const Record *dbRec);


};

#endif /* NewChromSweep_H */
