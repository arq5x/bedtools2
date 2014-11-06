/*****************************************************************************
  NewChromsweep.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "NewChromsweep.h"
#include "ContextIntersect.h"
#include "FileRecordMgr.h"

NewChromSweep::NewChromSweep(ContextIntersect *context)
:	_context(context),
 	_queryFRM(NULL),
 	_numDBs(_context->getNumDatabaseFiles()),
 	_queryRecordsTotalLength(0),
 	_databaseRecordsTotalLength(0),
 	_wasInitialized(false),
 	_currQueryRec(NULL),
 	_runToQueryEnd(false)
{
}


bool NewChromSweep::init() {
    
	//Create new FileRecordMgrs for the input files.
	//Open them, and get the first record from each.
	//otherwise, return true.
    _queryFRM = _context->getFile(_context->getQueryFileIdx());
    
    _dbFRMs.resize(_numDBs, NULL);
    for (int i=0; i < _numDBs; i++) {
    	_dbFRMs[i] = _context->getDatabaseFile(i);
    }

    _currDbRecs.resize(_numDBs, NULL);
    for (int i=0; i < _numDBs; i++) {
    	nextRecord(false, i);
    }

    _caches.resize(_numDBs);

    //determine whether to stop when the database end is hit, or keep going until the
    //end of the query file is hit as well.

    if (_context->getNoHit() || _context->getWriteCount() || _context->getWriteOverlap() || _context->getWriteAllOverlap() || _context->getLeftJoin()) {
    	_runToQueryEnd = true;
    }
    _wasInitialized = true;
    return true;
 }

void NewChromSweep::closeOut() {
	while (!_queryFRM->eof()) {
		nextRecord(true);
	}

    for (int i=0; i < _numDBs; i++) {
    	while (!_dbFRMs[i]->eof()) {
    		nextRecord(false, i);
    	}
    }
}

NewChromSweep::~NewChromSweep(void) {
	if (!_wasInitialized) {
		return;
	}
	_queryFRM->deleteRecord(_currQueryRec);
	_currQueryRec = NULL;


    for (int i=0; i < _numDBs; i++) {
    	_dbFRMs[i]->deleteRecord(_currDbRecs[i]);
    	_currDbRecs[i] = NULL;
    }

	_queryFRM->close();

   for (int i=0; i < _numDBs; i++) {
	   _dbFRMs[i]->close();
   }
}


void NewChromSweep::scanCache(int dbIdx, RecordKeyVector &retList) {
	recListIterType cacheIter = _caches[dbIdx].begin();
    while (cacheIter != _caches[dbIdx].end())
    {
    	const Record *cacheRec = cacheIter->value();
        if (_currQueryRec->sameChrom(cacheRec) && !(_currQueryRec->after(cacheRec))) {
            if (intersects(_currQueryRec, cacheRec)) {
                retList.push_back(cacheRec);
            } else break; // cacheRec is after the query rec, stop scanning.
            cacheIter = _caches[dbIdx].next();
        }
        else {
            cacheIter = _caches[dbIdx].deleteCurrent();
    		_dbFRMs[dbIdx]->deleteRecord(cacheRec);
        }
    }
}

void NewChromSweep::clearCache(int dbIdx)
{
	//delete all objects pointed to by cache
	recListType &cache = _caches[dbIdx];
	for (recListIterType iter = cache.begin(); iter != cache.end(); iter = cache.next()) {
		_dbFRMs[dbIdx]->deleteRecord(iter->value());
	}
	cache.clear();
}

void NewChromSweep::masterScan(RecordKeyVector &retList) {
	for (int i=0; i < _numDBs; i++) {
		if (dbFinished(i) || chromChange(i, retList)) {
			continue;
		} else {

			// scan the database cache for hits
			scanCache(i, retList);
			//skip if we hit the end of the DB
			// advance the db until we are ahead of the query. update hits and cache as necessary
			while (_currDbRecs[i] != NULL &&
					_currQueryRec->sameChrom(_currDbRecs[i]) &&
					!(_currDbRecs[i]->after(_currQueryRec))) {
				if (intersects(_currQueryRec, _currDbRecs[i])) {
					retList.push_back(_currDbRecs[i]);
				}
				if (_currQueryRec->after(_currDbRecs[i])) {
					_dbFRMs[i]->deleteRecord(_currDbRecs[i]);
					_currDbRecs[i] = NULL;
				} else {
					_caches[i].push_back(_currDbRecs[i]);
					_currDbRecs[i] = NULL;
				}
				nextRecord(false, i);
			}
		}
	}
}

bool NewChromSweep::chromChange(int dbIdx, RecordKeyVector &retList)
{
    // the files are on the same chrom
	if (_currDbRecs[dbIdx] == NULL || _currQueryRec->sameChrom(_currDbRecs[dbIdx])) {
		return false;
	}
	// the query is ahead of the database. fast-forward the database to catch-up.
	if (_currQueryRec->chromAfter(_currDbRecs[dbIdx])) {

		while (_currDbRecs[dbIdx] != NULL &&
				_currQueryRec->chromAfter(_currDbRecs[dbIdx])) {
			_dbFRMs[dbIdx]->deleteRecord(_currDbRecs[dbIdx]);
			nextRecord(false, dbIdx);
		}
		clearCache(dbIdx);
        return false;
    }
    // the database is ahead of the query.
    else {
        // 1. scan the cache for remaining hits on the query's current chrom.
		scanCache(dbIdx, retList);

        return true;
    }

	//control can't reach here, but compiler still wants a return statement.
	return true;
}


bool NewChromSweep::next(RecordKeyVector &retList) {
	retList.clearVector();
	if (_currQueryRec != NULL) {
		_queryFRM->deleteRecord(_currQueryRec);
	}

	if (!nextRecord(true)) return false; // query EOF hit


	if (allCurrDBrecsNull() && allCachesEmpty() && !_runToQueryEnd) {
		return false;
	}
	_currChromName = _currQueryRec->getChrName();

	masterScan(retList);

	if (_context->getSortOutput()) {
		retList.sortVector();
	}

	retList.setKey(_currQueryRec);
	return true;
}

bool NewChromSweep::nextRecord(bool query, int dbIdx) {
	if (query) {
		_currQueryRec = _queryFRM->getNextRecord();
		if (_currQueryRec != NULL) {
			_queryRecordsTotalLength += (unsigned long)(_currQueryRec->getEndPos() - _currQueryRec->getStartPos());
			return true;
		}
		return false;
	} else { //database
		Record *rec = _dbFRMs[dbIdx]->getNextRecord();
		_currDbRecs[dbIdx] = rec;
		if (rec != NULL) {
			_databaseRecordsTotalLength += (unsigned long)(rec->getEndPos() - rec->getStartPos());
			return true;
		}
		return false;
	}
}

bool NewChromSweep::intersects(const Record *rec1, const Record *rec2) const
{
	return rec1->sameChromIntersects(rec2, _context->getSameStrand(), _context->getDiffStrand(),
			_context->getOverlapFraction(), _context->getReciprocal());
}


bool NewChromSweep::allCachesEmpty() {
	for (int i=0; i < _numDBs; i++) {
		if (!_caches[i].empty()) {
			return false;
		}
	}
	return true;
}

bool NewChromSweep::allCurrDBrecsNull() {
	for (int i=0; i < _numDBs; i++) {
		if (_currDbRecs[i] != NULL) {
			return false;
		}
	}
	return true;
}

bool NewChromSweep::dbFinished(int dbIdx) {
	if (_currDbRecs[dbIdx] == NULL && _caches[dbIdx].empty()) {
		return true;
	}
	return false;
}
