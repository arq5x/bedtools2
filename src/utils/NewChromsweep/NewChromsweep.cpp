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

NewChromSweep::NewChromSweep(ContextIntersect *context,
                       bool useMergedIntervals)
:	_context(context),
 	_queryFRM(NULL),
 	_databaseFRM(NULL),
 	_useMergedIntervals(false),
 	_queryRecordsTotalLength(0),
 	_databaseRecordsTotalLength(0),
 	_wasInitialized(false),
 	_currQueryRec(NULL),
 	_currDatabaseRec(NULL),
 	_runToQueryEnd(false)

{
}


bool NewChromSweep::init() {
    
	//Create new FileRecordMgrs for the two input files.
	//Open them, and get the first record from each.
	//if any of that goes wrong, return false;
	//otherwise, return true.
    _queryFRM = _context->getFile(_context->getQueryFileIdx());
    _databaseFRM = _context->getFile(_context->getDatabaseFileIdx());
    
    nextRecord(false);
//    if (_currDatabaseRec == NULL) {
//    	return false;
//    }

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
	while (!_databaseFRM->eof()) {
		nextRecord(false);
	}
}

NewChromSweep::~NewChromSweep(void) {
	if (!_wasInitialized) {
		return;
	}
	_queryFRM->deleteRecord(_currQueryRec);
	_currQueryRec = NULL;

	_databaseFRM->deleteRecord(_currDatabaseRec);
	_currDatabaseRec = NULL;

	_queryFRM->close();
	_databaseFRM->close();
}


void NewChromSweep::scanCache() {
	recListIterType cacheIter = _cache.begin();
    while (cacheIter != _cache.end())
    {
    	const Record *cacheRec = cacheIter->value();
        if (_currQueryRec->sameChrom(cacheRec) && !(_currQueryRec->after(cacheRec))) {
            if (intersects(_currQueryRec, cacheRec)) {
                _hits.push_back(cacheRec);
            }
            cacheIter = _cache.next();
        }
        else {
            cacheIter = _cache.deleteCurrent();
    		_databaseFRM->deleteRecord(cacheRec);
        }
    }
}

void NewChromSweep::clearCache()
{
	//delete all objects pointed to by cache
	for (recListIterType iter = _cache.begin(); iter != _cache.end(); iter = _cache.next()) {
		_databaseFRM->deleteRecord(iter->value());
	}
	_cache.clear();
}

bool NewChromSweep::chromChange()
{
    // the files are on the same chrom
	if (_currDatabaseRec == NULL || _currQueryRec->sameChrom(_currDatabaseRec)) {
		return false;
	}
	// the query is ahead of the database. fast-forward the database to catch-up.
	if (_currQueryRec->chromAfter(_currDatabaseRec)) {

		while (_currDatabaseRec != NULL &&
				_currQueryRec->chromAfter(_currDatabaseRec)) {
			_databaseFRM->deleteRecord(_currDatabaseRec);
			nextRecord(false);
		}
		clearCache();
        return false;
    }
    // the database is ahead of the query.
    else {
        // 1. scan the cache for remaining hits on the query's current chrom.
        if (_currQueryRec->getChrName() == _currChromName)
        {
            scanCache();
        }
        // 2. fast-forward until we catch up and report 0 hits until we do.
        else if (_currQueryRec->chromBefore(_currDatabaseRec))
        {
            clearCache();
        }

        return true;
    }
}


bool NewChromSweep::next(RecordKeyList &next) {
	if (_currQueryRec != NULL) {
		_queryFRM->deleteRecord(_currQueryRec);
	}
	nextRecord(true);
	if (_currQueryRec == NULL) { //eof hit!
		return false;
	}

	if (_currDatabaseRec == NULL && _cache.empty() && !_runToQueryEnd) {
		return false;
	}
	_hits.clear();
	_currChromName = _currQueryRec->getChrName();
	// have we changed chromosomes?
	if (!chromChange()) {
		// scan the database cache for hits
		scanCache();
		//skip if we hit the end of the DB
		// advance the db until we are ahead of the query. update hits and cache as necessary
		while (_currDatabaseRec != NULL &&
				_currQueryRec->sameChrom(_currDatabaseRec) &&
				!(_currDatabaseRec->after(_currQueryRec))) {
			if (intersects(_currQueryRec, _currDatabaseRec)) {
				_hits.push_back(_currDatabaseRec);
			}
			if (_currQueryRec->after(_currDatabaseRec)) {
				_databaseFRM->deleteRecord(_currDatabaseRec);
				_currDatabaseRec = NULL;
			} else {
				_cache.push_back(_currDatabaseRec);
				_currDatabaseRec = NULL;
			}
			nextRecord(false);
		}
	}
	next.setKey(_currQueryRec);
	next.setListNoCopy(_hits);
	return true;
}

void NewChromSweep::nextRecord(bool query) {
	if (query) {
//		if (!_context->getUseMergedIntervals()) {
			_currQueryRec = _queryFRM->getNextRecord();
//		} else {
//			_currQueryRec = _queryFRM->allocateAndGetNextMergedRecord(_context->getSameStrand() ? FileRecordMgr::SAME_STRAND_EITHER : FileRecordMgr::ANY_STRAND);
//		}
		if (_currQueryRec != NULL) {
			_queryRecordsTotalLength += (unsigned long)(_currQueryRec->getEndPos() - _currQueryRec->getStartPos());
		}
	} else { //database
//		if (!_context->getUseMergedIntervals()) {
			_currDatabaseRec = _databaseFRM->getNextRecord();
//		} else {
//			_currDatabaseRec = _databaseFRM->allocateAndGetNextMergedRecord(_context->getSameStrand() ? FileRecordMgr::SAME_STRAND_EITHER : FileRecordMgr::ANY_STRAND);
//		}
		if (_currDatabaseRec != NULL) {
			_databaseRecordsTotalLength += (unsigned long)(_currDatabaseRec->getEndPos() - _currDatabaseRec->getStartPos());
		}
	}
}

bool NewChromSweep::intersects(const Record *rec1, const Record *rec2) const
{
	return rec1->sameChromIntersects(rec2, _context->getSameStrand(), _context->getDiffStrand(),
			_context->getOverlapFraction(), _context->getReciprocal());
}
