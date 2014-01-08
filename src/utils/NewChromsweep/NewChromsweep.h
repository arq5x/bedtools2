/*****************************************************************************
  NewChromsweepBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef NEW_CHROMSWEEP_H
#define NEW_CHROMSWEEP_H

using namespace std;

#include <string>
#include "BTlist.h"
#include "RecordKeyList.h"
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
    
    typedef BTlist<const Record *> recListType;
    typedef const BTlistNode<const Record *> * recListIterType;
    // loads next (a pair) with the current query and it's overlaps
    bool next(RecordKeyList &next);
    
    //MUST call this method after sweep if you want the
    //getTotalRecordLength methods to return the whole length of the
    //record files, rather than just what was used by sweep.
    void closeOut();

    unsigned long getQueryTotalRecordLength() { return _queryRecordsTotalLength; }
    unsigned long getDatabaseTotalRecordLength() { return _databaseRecordsTotalLength; }

    // Usage:
    //     NewChromSweep sweep = NewChromSweep(queryFileName, databaseFileName);
    //     RecordKeyList hits;
    //     while (sweep.next(hits))
    //     {
    //        processHits(hits);
    //     }
    // 	   getQueryTotalRecordLength()
    //     getDatabaseTotalRecordLength()

private:
    ContextIntersect *_context;
    FileRecordMgr *_queryFRM;
    FileRecordMgr *_databaseFRM;

    bool _useMergedIntervals;

    unsigned long _queryRecordsTotalLength;
    unsigned long _databaseRecordsTotalLength;

    bool _wasInitialized;

    // a cache of still active features from the database file
    recListType _cache;
    // the set of hits in the database for the current query
    recListType _hits;

    // the current query and db features.
    Record * _currQueryRec;
    Record *_currDatabaseRec;
    // a cache of the current chrom from the query. used to handle chrom changes.
    QuickString _currChromName;
    bool _runToQueryEnd;

    void nextRecord(bool query); //true fetches next query record, false fetches next db record.
    void nextDatabase();
    
    void scanCache();
    void clearCache();
    bool chromChange();
    bool intersects(const Record *rec1, const Record *rec2) const;

};

#endif /* NewChromSweep_H */
