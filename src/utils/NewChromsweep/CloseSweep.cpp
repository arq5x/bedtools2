/*
 * CloseSweep.cpp
 *
 *  Created on: Sep 25, 2014
 *      Author: nek3d
 */

#include "CloseSweep.h"
#include "ContextClosest.h"

RecDistList::RecDistList(int maxSize)
:  _kVal(maxSize),
   _empty(true),
   _currNumIdxs(0),
   _totalRecs(0)

    {

	_allRecs.resize(_kVal);
	for (int i=0; i < _kVal; i++) {
		_allRecs[i] = new elemsType();
	}
	_distIndex = new indexType[_kVal];

	clear();
}

RecDistList::~RecDistList() {
	for (int i=0; i < _kVal; i++) {
		delete _allRecs[i];
	}
	delete _distIndex;
	_distIndex = NULL;

}

void RecDistList::clear() {
	for (int i=0; i < _kVal; i++) {
		_allRecs[i]->clear();
		_allRecs[i]->reserve(16);
	}
	for (int i=0; i < _kVal; i++) {
		_distIndex[i].first = -1;
		_distIndex[i].second = -1;
	}
	_currNumIdxs = 0;
	_empty = true;
	_totalRecs = 0;
}

bool RecDistList::addRec(CHRPOS dist, Record *record, chromDirType chromDir) {
	CHRPOS newPos = 0;
	bool mustAppend = false;
	int useElemIdx = 0;

	if (dist > getMaxDist()) {
		//dist is bigger than any currently contained.
		if (_currNumIdxs == _kVal) {
			return false;
		}
		mustAppend = true;
		useElemIdx = _currNumIdxs;
		newPos = _currNumIdxs;
	}
	if (find(dist, newPos)) {
		_allRecs[_distIndex[newPos].second]->push_back(elemPairType(chromDir, record));
		_totalRecs++;
		return true;
	}
	if (!mustAppend) {
		//smaller than maxDist, and doesn't currently exist. must insert
		//newPos now is the insertion point.
		CHRPOS startShiftPos = 0;
		if (_currNumIdxs == _kVal) {
			//already full. must remove oldest max.
			//determine which vector it was using
			//so we can re-use it.
			startShiftPos = _kVal-1;
			useElemIdx = _distIndex[startShiftPos].second;
		} else {
			//can add a new element
			startShiftPos = _currNumIdxs;
			useElemIdx = _currNumIdxs;
			_currNumIdxs++;
		}
		for (CHRPOS i=startShiftPos; i > newPos; i--) {
			_distIndex[i].first = _distIndex[i-1].first;
			_distIndex[i].second = _distIndex[i-1].second;
		}
	} else {
		_currNumIdxs++;
	}
	_allRecs[useElemIdx]->clear();
	_allRecs[useElemIdx]->reserve(16);
	_allRecs[useElemIdx]->push_back(elemPairType(chromDir, record));

	_distIndex[newPos].first = (int)dist;
	_distIndex[newPos].second = useElemIdx;
	_empty = false;
	_totalRecs++;
	return true;
}

//if true, pos will be the idx the distance is at.
//if false, pos will be the idx to insert at.
bool RecDistList::find(CHRPOS dist, CHRPOS &pos) const {
	CHRPOS lbound=0, ubound=_currNumIdxs-1, currVal =0;
	pos = 0;
	while(lbound <= ubound)
	{
		pos = (lbound + ubound) / 2;
		currVal = _distIndex[pos].first;
		if (currVal == dist) {
			return true;
		}
		if (dist > currVal) {
			lbound = pos + 1;
		} else {
			ubound = pos -1;
		}
	}
	pos = (ubound == -1 ) ? 0 : (lbound == _currNumIdxs ? _currNumIdxs : lbound);
	return false;
}

CHRPOS RecDistList::getMaxLeftEndPos() const {

	if (_empty) return -1;

	CHRPOS maxDist =_distIndex[_currNumIdxs-1].first;

  const elemsType *elems = _allRecs[_distIndex[_currNumIdxs-1].second];
	for (int i=0; i < (int)elems->size(); i++) {
    const elemPairType & elem = (*elems)[i];
		if (elem.first == LEFT) {
			return maxDist;
		}
	}
	return -1;

}


CloseSweep::CloseSweep(ContextClosest *context)
:	NewChromSweep(context),
 	_context(context),
 	_kClosest(_context->getNumClosestHitsWanted()),
	_sameStrand(false),
	_diffStrand(false),

	_refDist(false),
	_aDist(false),
	_bDist(false),

	_ignoreUpstream(false),
	_ignoreDownstream(false),

	_qForward(false),
	_qReverse(false),
	_dbForward(false),
	_dbReverse(false),

	_tieMode(ContextClosest::ALL_TIES),
	_firstTie(false),
	_lastTie(false),
	_allTies(false)
 	{

	_minUpstreamRecs.resize(_numDBs, NULL);
	_minDownstreamRecs.resize(_numDBs, NULL);
	_overlapRecs.resize(_numDBs, NULL);
	_maxPrevLeftClosestEndPos.resize(_numDBs, 0);
	_maxPrevLeftClosestEndPosReverse.resize(_numDBs, 0);

	for (int i=0; i < _numDBs; i++) {
		_minUpstreamRecs[i] = new RecDistList(_kClosest);
		_minDownstreamRecs[i] = new RecDistList(_kClosest);
		_overlapRecs[i] = new RecDistList(_kClosest);
	}
}

CloseSweep::~CloseSweep(void) {
	for (int i=0; i < _numDBs; i++) {
		delete _minUpstreamRecs[i];
		delete _minDownstreamRecs[i];
		delete _overlapRecs[i];
	}
}

bool CloseSweep::init() {

    bool retVal =  NewChromSweep::init();
    _runToQueryEnd = true;

	// Some abbreviations to make the code less miserable.
	_sameStrand = _context->getSameStrand();
	_diffStrand = _context->getDiffStrand();

	_refDist = _context->getStrandedDistMode() == ContextClosest::REF_DIST;
	_aDist = _context->getStrandedDistMode() == ContextClosest::A_DIST;
	_bDist = _context->getStrandedDistMode() == ContextClosest::B_DIST;

	_ignoreUpstream = _context->ignoreUpstream();
	_ignoreDownstream = _context->ignoreDownstream();

	_tieMode = _context->getTieMode();
	_firstTie = _tieMode == ContextClosest::FIRST_TIE;
	_lastTie = _tieMode == ContextClosest::LAST_TIE;
	_allTies = _tieMode == ContextClosest::ALL_TIES;


    return retVal;
 }

void CloseSweep::masterScan(RecordKeyVector &retList) {

	_qForward = _currQueryRec->getStrandVal() == Record::FORWARD;
	_qReverse = _currQueryRec->getStrandVal() == Record::REVERSE;

	if (_currQueryChromName != _prevQueryChromName) testChromOrder(_currQueryRec);
	if (_context->reportDistance()) {
		_finalDistances.clear();
	}

	for (int i=0; i < _numDBs; i++) {

		//first clear out everything from the previous scan
		_minUpstreamRecs[i]->clear();
		_minDownstreamRecs[i]->clear();
		_overlapRecs[i]->clear();

		if (dbFinished(i) || chromChange(i, retList, true)) {
			continue;
		} else {

			// scan the database cache for hits
			scanCache(i, retList);

			// skip if we hit the end of the DB
			// advance the db until we are ahead of the query. update hits and cache as necessary
			bool stopScanning = false;
			while (_currDbRecs[i] != NULL &&
					_currQueryRec->sameChrom(_currDbRecs[i]) &&
					!stopScanning) {
				if (considerRecord(_currDbRecs[i], i, stopScanning) == DELETE) {
					_dbFRMs[i]->deleteRecord(_currDbRecs[i]);
					_currDbRecs[i] = NULL;
				} else {
					_caches[i].push_back(_currDbRecs[i]);
					_currDbRecs[i] = NULL;
				}
				nextRecord(false, i);
			}
		}
		finalizeSelections(i, retList);
	}
	checkMultiDbs(retList);
}

void CloseSweep::scanCache(int dbIdx, RecordKeyVector &retList) {
	recListIterType cacheIter = _caches[dbIdx].begin();
    while (cacheIter != _caches[dbIdx].end())
    {
        Record *cacheRec = cacheIter->value();
    	bool stopScanning = false;
    	if (considerRecord(cacheRec, dbIdx, stopScanning) == DELETE) {
            cacheIter = _caches[dbIdx].deleteCurrent();
    		_dbFRMs[dbIdx]->deleteRecord(cacheRec);
    	} else {
            cacheIter = _caches[dbIdx].next();
    	}
    	if (stopScanning) break;
    }
}


CloseSweep::rateOvlpType CloseSweep::considerRecord(Record *cacheRec, int dbIdx, bool &stopScanning) {

	// Determine whether the hit and query intersect, and if so, what to do about it.
	_dbForward = cacheRec->getStrandVal() == Record::FORWARD;
	_dbReverse = cacheRec->getStrandVal() == Record::REVERSE;
	CHRPOS currDist = 0;

	if (intersects(_currQueryRec, cacheRec)) {

		// HIT INTERSECTS QUERY
		return tryToAddRecord(cacheRec, 0, dbIdx, stopScanning, OVERLAP, INTERSECT);

	} else if (cacheRec->after(_currQueryRec)) {

		// HIT IS TO THE RIGHT OF THE QUERY.

		 currDist = (cacheRec->getStartPos() - _currQueryRec->getEndPos()) + 1;
		 if (_context->signDistance()) {
			 if ((_aDist && _qReverse) ||
				 (_bDist && _dbForward))
				 {
					 // hit is "upstream" of A
					 return tryToAddRecord(cacheRec, abs(currDist), dbIdx, stopScanning, RIGHT, UPSTREAM);
				 }
		 }
		 // HIT IS DOWNSTREAM.
		 return tryToAddRecord(cacheRec, abs(currDist), dbIdx, stopScanning, RIGHT, DOWNSTREAM);
	 } else if (_currQueryRec->after(cacheRec)){

		 // HIT IS TO THE LEFT OF THE QUERY.
		 currDist = (_currQueryRec->getStartPos() - cacheRec->getEndPos()) + 1;
		 if (_context->signDistance()) {
			 if ((_aDist && _qReverse) ||
				 (_bDist && _dbForward))
			 {
				 // HIT IS DOWNSTREAM.
				 return tryToAddRecord(cacheRec, abs(currDist), dbIdx, stopScanning, LEFT, DOWNSTREAM);
			 }
		 }
		 // hit is "UPSTREAM" of A
		 return tryToAddRecord(cacheRec, abs(currDist), dbIdx, stopScanning, LEFT, UPSTREAM);
	 }
	return IGNORE;
}

void CloseSweep::finalizeSelections(int dbIdx, RecordKeyVector &retList) {

	// Determine whether the overlaps, upstream, or downstream records have
	// the k closest hits.

	// The first thing to do is set the leftmost end pos used by records on the left.
	// This will control when the cache is purged during the next query's sweep.
	setLeftClosestEndPos(dbIdx);



	RecDistList *upRecs = _minUpstreamRecs[dbIdx];
	RecDistList *downRecs = _minDownstreamRecs[dbIdx];
	RecDistList *overlaps = _overlapRecs[dbIdx];
	RecDistList::constIterType upIter = upRecs->begin();
	RecDistList::constIterType downIter = downRecs->begin();
	CHRPOS upDist = INT_MAX;
	CHRPOS downDist = INT_MAX;

	int totalHitsUsed = 0;


	// If forcing upstream, use those first.
	if (_context->forceUpstream()) {
		//add upstream recs until all are used or K hits taken.
		while (upIter != upRecs->end() && totalHitsUsed < _kClosest) {
			upDist = upRecs->currDist(upIter);
			totalHitsUsed += addRecsToRetList(upRecs->allElems(upIter), 0 - upDist, retList);
			upIter++;
		}
	}

	// If forcing downstream, use those first/next.
	if (_context->forceDownstream()) {
		while (downIter != downRecs->end() && totalHitsUsed < _kClosest) {
			downDist = downRecs->currDist(downIter);
			totalHitsUsed += addRecsToRetList(downRecs->allElems(downIter), downDist, retList);
			downIter++;
		}
	}


	//start with the overlaps, which will all have distance zero.
	if (totalHitsUsed < _kClosest && !overlaps->empty()) {
		//there are overlaps.
		totalHitsUsed += addRecsToRetList(overlaps->allElems(overlaps->begin()), 0, retList);
	}

	//now check the upstream and downstream recs, starting with the beginning of each,
	//and grabbing whichever has the closest records as we work our way out to the
	//max dist. Continue until we have K records, or the reclists are empty.

	while (totalHitsUsed < _kClosest) {
		upDist = INT_MAX;
		downDist = INT_MAX;

		if (upIter != upRecs->end()) {
			upDist = upRecs->currDist(upIter);
		}
		if (downIter != downRecs->end()) {
			downDist = downRecs->currDist(downIter);
		}

		//stop if no hits are left to consider
		if ((upDist == INT_MAX) && (downDist == INT_MAX)) break;

		bool tie = upDist == downDist;
		bool usedUp = false;
		bool usedDown = false;
		if (upDist < downDist || (tie && !_lastTie)) {
			totalHitsUsed += addRecsToRetList(upRecs->allElems(upIter), 0 - upDist, retList);
			upIter++;
			usedUp = true;
		}
		if (downDist < upDist || (tie && !_firstTie)) {
			totalHitsUsed += addRecsToRetList(downRecs->allElems(downIter), downDist, retList);
			downIter++;
			usedDown = true;
		}
		if (tie) {
			// If there was a tie, but we didn't use both elements because of the tie mode,
			// then we still have to increment the iterator of the unused element so it
			// isn't used later.
			if (usedUp && !usedDown) {
				downIter++;
			} else if (usedDown && !usedUp) {
				upIter++;
			}
		}
	}

}



int CloseSweep::addRecsToRetList(RecDistList::elemsType *recs, CHRPOS currDist, RecordKeyVector &retList) {

	int hitsUsed = 0;
	int numRecs = (int)recs->size(); //just to clean the code some.


	if (_firstTie) {
		addSingleRec(recs->at(0).second, currDist, hitsUsed, retList);
		return 1;

	} else if (_lastTie) {
		addSingleRec(recs->at(numRecs-1).second, currDist, hitsUsed, retList);
		return 1;

	} else { //tieMode == ALL_TIES

		for (int i=0; i < numRecs; i++) {
			addSingleRec(recs->at(i).second, currDist, hitsUsed, retList);
		}
		return numRecs;
	}
}

void CloseSweep::addSingleRec(Record *rec, CHRPOS currDist, int &hitsUsed, RecordKeyVector &retList) {
	retList.push_back(rec);
	_finalDistances.push_back(currDist);
	hitsUsed++;
}

void CloseSweep::checkMultiDbs(RecordKeyVector &retList) {
//	//can skip this method if there's only one DB, or if we are
//	//resolving closest hits for each db instead of all of them
	if (_context->getMultiDbMode() != ContextClosest::ALL_DBS ||  _numDBs == 1) return;


	// Get the K closest hits among multiple databases,
	// while not counting ties more than once if the tieMode
	// is "first" or "last".
	// Start by entering  all hits and their absolute distances
	// into a vector of distance tuples, then sort it.

	vector<distanceTuple> copyDists;
	int numHits = (int)retList.size();
	copyDists.resize(numHits);
	CHRPOS i=0;
	for (RecordKeyVector::iterator_type iter = retList.begin(); iter != retList.end(); iter++) {
		CHRPOS dist = _finalDistances[i];
		copyDists[i]._dist = abs(dist);
		copyDists[i]._rec = *iter;
		copyDists[i]._isNeg = dist < 0;
		i++;
	}

	// sort the hits by distance
	sort(copyDists.begin(), copyDists.end(), DistanceTupleSortAscFunctor());

	//now we want to build a map telling us what distances are tied,
	//and how many of each of these there are. Use a map<int, int>,
	//where the key is a distance (in absolute value) and the value
	//is the number of ties that that distance has.
	map<CHRPOS, CHRPOS> ties;
	for (vector<distanceTuple>::iterator i = copyDists.begin(); i != copyDists.end(); ++i)
    	++ties[i->_dist];

	// Clear the original list and distances, and re-populate
	// until we have the desired number of hits, skipping
	// over any unwanted ties.
	retList.clearVector();
	_finalDistances.clear();

	int hitsUsed = 0;
	for (i=0; i < numHits && hitsUsed < _kClosest; i++) {
		CHRPOS dist = copyDists[i]._dist;
		bool isNeg = copyDists[i]._isNeg;
		//see if this distance is tied with any other
		map<CHRPOS, CHRPOS>::iterator iter = ties.find(dist);
		if (iter != ties.end()) {
			//tie was found
			CHRPOS numTies = iter->second;
			if (!_allTies) {
				if (_firstTie) {
					//just add the first of the ties
					addSingleRec(copyDists[i]._rec, (isNeg ? 0 - dist : dist), hitsUsed, retList);
					i += numTies - 1; // use first, then skip ahead by the number of ties, minus 1 because
					//loop is about to be incremented
				} else { //tieMode == LAST_TIE. Just add the last of the ties.
					i += numTies -1;
					dist = copyDists[i]._dist;
					isNeg = copyDists[i]._isNeg;
					addSingleRec(copyDists[i]._rec, (isNeg ? 0 - dist : dist), hitsUsed, retList);
				}
			} else {
				// tieMode is ALL_TIES, use all hits.
				for (int j = i; j < i + numTies; j++) {
					dist = copyDists[j]._dist;
					isNeg = copyDists[j]._isNeg;
					addSingleRec(copyDists[j]._rec, (isNeg ? 0 - dist : dist), hitsUsed, retList);
				}
				i += numTies - 1; //skip ahead by the number of ties, minus 1 because
				//loop is about to be incremented
			}
		} else {
			addSingleRec(copyDists[i]._rec, (isNeg ? 0 - dist : dist), hitsUsed, retList);
		}
	}
}

bool CloseSweep::chromChange(int dbIdx, RecordKeyVector &retList, bool wantScan)
{
	Record *dbRec = _currDbRecs[dbIdx];

	bool haveQuery = _currQueryRec != NULL;
	bool haveDB = dbRec != NULL;

	if (haveQuery && _currQueryChromName != _prevQueryChromName) {
		_context->testNameConventions(_currQueryRec);
		testChromOrder(_currQueryRec);
	}

	if (haveDB) {
		_context->testNameConventions(dbRec);
		testChromOrder(dbRec);
	}

    // the files are on the same chrom
	if (haveQuery && (!haveDB || _currQueryRec->sameChrom(dbRec))) {

		//if this is the first time the query's chrom is ahead of the chrom that was in this cache,
		//then we have to clear the cache.
		if (!_caches[dbIdx].empty() && queryChromAfterDbRec(_caches[dbIdx].begin()->value())) {
			clearCache(dbIdx);
			clearClosestEndPos(dbIdx);
		}
		return false;
	}

	if (!haveQuery || !haveDB) return false;

	if (!_caches[dbIdx].empty() && (_caches[dbIdx].begin()->value()->sameChrom(_currQueryRec))) {
		//the newest DB record's chrom is ahead of the query, but the cache still
		//has old records on that query's chrom
		scanCache(dbIdx, retList);
		finalizeSelections(dbIdx, retList);
		return true;
	}


	// the query is ahead of the database. fast-forward the database to catch-up.
	if (queryChromAfterDbRec(dbRec)) {
		string oldDbChrom(dbRec->getChrName());

		while (dbRec != NULL &&
				queryChromAfterDbRec(dbRec)) {
			_dbFRMs[dbIdx]->deleteRecord(dbRec);
			if (!nextRecord(false, dbIdx)) break;
			dbRec =  _currDbRecs[dbIdx];
			const string &newDbChrom = dbRec->getChrName();
			if (newDbChrom != oldDbChrom) {
				testChromOrder(dbRec);
				oldDbChrom = newDbChrom;
			}
		}
		clearCache(dbIdx);
		clearClosestEndPos(dbIdx);
        return false;
    }
    // the database is ahead of the query.
    else {
        // 1. scan the cache for remaining hits on the query's current chrom.
		if (wantScan) scanCache(dbIdx, retList);

        return true;
    }

	//control can't reach here, but compiler still wants a return statement.
	return true;
}


void CloseSweep::setLeftClosestEndPos(int dbIdx)
{

  //try to determine max end pos of hits to left of query.
  //first check for exceptions due to options that cause
  //complete rejection of all left side hits.

  //no records found to left of query,
  //can't set purge point, except for some special cases.
  purgeDirectionType purgeDir = purgePointException();
  if (purgeDir == BOTH || purgeDir == FORWARD_ONLY) {
    _maxPrevLeftClosestEndPos[dbIdx] = max(_currQueryRec->getStartPos(), _maxPrevLeftClosestEndPos[dbIdx]);
  }
  if (purgeDir == BOTH || purgeDir == REVERSE_ONLY) {
    _maxPrevLeftClosestEndPosReverse[dbIdx] = max(_currQueryRec->getStartPos(), _maxPrevLeftClosestEndPosReverse[dbIdx]);
  }
  if (purgeDir != NEITHER) return;

  RecDistList *upRecs = _minUpstreamRecs[dbIdx];
  RecDistList *downRecs = _minDownstreamRecs[dbIdx];

  CHRPOS upDist = upRecs->getMaxLeftEndPos();
  CHRPOS downDist = downRecs->getMaxLeftEndPos();


  if (upDist == -1 && downDist == -1) return;

	CHRPOS leftMostEndPos = (_currQueryRec->getStartPos() - max(upDist, downDist)) +1;
	if ((!_sameStrand && !_diffStrand) ||
		(_sameStrand && _qForward) ||
		(_diffStrand && _qReverse))  {

		_maxPrevLeftClosestEndPos[dbIdx] = max(leftMostEndPos, _maxPrevLeftClosestEndPos[dbIdx]);
	} else {
		_maxPrevLeftClosestEndPosReverse[dbIdx] = max(leftMostEndPos, _maxPrevLeftClosestEndPosReverse[dbIdx]);
	}
}

bool CloseSweep::beforeLeftClosestEndPos(int dbIdx, Record *rec)
{
	CHRPOS recEndPos = rec->getEndPos();
	CHRPOS prevPos = _maxPrevLeftClosestEndPos[dbIdx];
	CHRPOS prevPosReverse = _maxPrevLeftClosestEndPosReverse[dbIdx];

	if (!_sameStrand && !_diffStrand) {
		return recEndPos < prevPos;
	} else {
		if (rec->getStrandVal() == Record::FORWARD) {
			return recEndPos < prevPos;
		} else {
			return recEndPos < prevPosReverse;
		}
  }
	return false;
}

void CloseSweep::clearClosestEndPos(int dbIdx)
{
	_maxPrevLeftClosestEndPos[dbIdx] = 0;
	_maxPrevLeftClosestEndPosReverse[dbIdx] = 0;
}

CloseSweep::rateOvlpType CloseSweep::tryToAddRecord(Record *cacheRec, CHRPOS dist, int dbIdx, bool &stopScanning, chromDirType chromDir, streamDirType streamDir) {

	//
	// Decide whether to ignore hit
	//
	// If want same strand, and they're unknown or not the same, ignore
	// If we want diff strand, and they're unknown or different, ignore.
	// If we want diff names and they're the same, ignore
	// If stream is unwanted, ignore.
	bool hasUnknownStrands = (_currQueryRec->getStrandVal() == Record::UNKNOWN || cacheRec->getStrandVal() == Record::UNKNOWN);
	bool wantedSameNotSame = (_sameStrand && (hasUnknownStrands || (_currQueryRec->getStrandVal() != cacheRec->getStrandVal())));
	bool wantedDiffNotDiff = (_diffStrand && (hasUnknownStrands || (_currQueryRec->getStrandVal() == cacheRec->getStrandVal())));
	bool badStrand = wantedSameNotSame || wantedDiffNotDiff;
	bool badNames = (_context->diffNames() && (cacheRec->getName() == _currQueryRec->getName()));
	bool badStream = (streamDir == UPSTREAM ? _ignoreUpstream : (streamDir == DOWNSTREAM ? _ignoreDownstream :  _context->ignoreOverlaps()));

	bool shouldIgnore = badStrand || badNames || badStream;


	// You would think ignoring it means we could stop here, but even then,
	// hits on the left may need to be purged from the cache,
	// and hits on the right tell us when to stop scanning the cache.


	//establish which set of hits we are looking at.
	RecDistList *useList = (streamDir == UPSTREAM ? _minUpstreamRecs[dbIdx] : (streamDir == INTERSECT ? _overlapRecs[dbIdx] : _minDownstreamRecs[dbIdx]));


	if (chromDir == OVERLAP && !shouldIgnore) {
		useList->addRec(0, cacheRec, RecDistList::OVERLAP);
	}

	else if (chromDir == LEFT) {
		if (beforeLeftClosestEndPos(dbIdx, cacheRec)) {
			return DELETE;
		}
		if (!shouldIgnore) {
			useList->addRec(dist, cacheRec, RecDistList::LEFT);
		}
	}

	else if (chromDir == RIGHT) {
		//hit is to the right of query. Need to know when to stop scanning.
		// if NOT ignoring:
		//		if we're able to add it to the useList, definitely DON'T stop,
		//		hit was valid, so next one could be too.
		//		if we're UNABLE to add it to the uselist, definitely DO stop,
		//		because useList is full.
		// but if we ARE ignoring:
		//      can only stop in the rare cases where hits to the right
		//		are ALWAYS getting ignored, regardless of query and cache strand
		//      Otherwise, always continue.
		if (!shouldIgnore) {
			if (!useList->addRec(dist, cacheRec, RecDistList::RIGHT)) {
				stopScanning = true;
			}
		} else { // hit is ignored
			if (allHitsRightOfQueryIgnored()) {
				stopScanning = true;
			}
		}
	}
	return IGNORE;
}

bool CloseSweep::allHitsRightOfQueryIgnored() {
	return ((_refDist && _ignoreDownstream) ||
			(_aDist && ((_ignoreUpstream && _qReverse) || (_ignoreDownstream && _qForward))) ||
			(_bDist && ((_ignoreDownstream && ((_qForward && _diffStrand) || (_qReverse && _sameStrand))) ||
					((_ignoreUpstream && ((_qReverse && _diffStrand) || (_qForward && _sameStrand)))))));
}


CloseSweep::purgeDirectionType CloseSweep::purgePointException() {

	// Normally, we can't set a cache purge point if there are no
	// records to the left of the query.

	// This method will detect use cases that cause all of the left side
	// hits to have been rejected, and tell us whether to purge the forward
	// cache, reverse cache, both, or neither.

	purgeDirectionType purgeDir = NEITHER;

	if (_ignoreUpstream && _ignoreDownstream) {
		purgeDir = BOTH;
	}
	else if (_ignoreUpstream) {
		if (_refDist) {
			purgeDir = BOTH;
		}
		else if (_aDist && _qForward) {
			if (_sameStrand) {
				purgeDir = FORWARD_ONLY;
			}
			else if (_diffStrand) {
				purgeDir = REVERSE_ONLY;
			}
		}
		else if (_bDist) {
			  if (_qForward && _diffStrand) purgeDir = REVERSE_ONLY;
			  else if (_qReverse && _sameStrand) purgeDir = REVERSE_ONLY;
		}
	}
	else if (_ignoreDownstream) {
		 // if refDist, do nothing. left hits can't be downstream.
		 if (_aDist) {
			//if qForward, do nothing. left hits can't be downstream.
			if (_qReverse) {
				if (_sameStrand) purgeDir = REVERSE_ONLY;
				else if (_diffStrand) purgeDir = FORWARD_ONLY;
			}
		} else if (_bDist) {
			if (_qForward && _sameStrand) purgeDir = FORWARD_ONLY;
			else if (_qReverse && _diffStrand) purgeDir = FORWARD_ONLY;
		}
	}

	return purgeDir;

}
