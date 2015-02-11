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
   _totalRecs(0) {

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

bool RecDistList::addRec(int dist, const Record *record, chromDirType chromDir) {
	int newPos = 0;
	bool mustAppend = false;
	int useElemIdx = 0;

	if (dist > getMaxDist()) {
		//dist is bigger than any currently contained.
		if (_currNumIdxs == _kVal) {
			//already full with smaller distances
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
		int startShiftPos = 0;
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
		for (int i=startShiftPos; i > newPos; i--) {
			_distIndex[i].first = _distIndex[i-1].first;
			_distIndex[i].second = _distIndex[i-1].second;
		}
	} else {
		_currNumIdxs++;
	}
	_allRecs[useElemIdx]->clear();
	_allRecs[useElemIdx]->reserve(16);
	_allRecs[useElemIdx]->push_back(elemPairType(chromDir, record));

	_distIndex[newPos].first = dist;
	_distIndex[newPos].second = useElemIdx;
	_empty = false;
	_totalRecs++;
	return true;
}

//bool RecDistList::addRec(int dist, const Record *record, chromDirType chromDir) {
//	//if unique size is >= maxRecNum, that means the collection is
//	//full. In that case, if the new rec has greater distance than
//	// any in the collection, ignore it. If the distance isn't
//	// greater, add it, removing the previously greater one
//	// if need be.
//	// If the collection isn't full, just the new.
//
//	if (uniqueSize() < _kVal || exists(dist)) {
//		insert(dist, record, chromDir);
//		return true;
//	} else {
//		//find previous greatest distance
//		revIterType rIter = _recs.rbegin();
//		int currMaxDist = rIter->first;
//		if (dist > currMaxDist) {
//			return false;
//		}
//		//Now we know dist is less than currMax.
//		//new dist is smaller. Erase old max.
//		delete rIter->second;
//
//		//apparently you can't erase a reverse iterator, so we'll
//		//make a normal one pointing at the end, then back it  up.
//		iterType delIter = _recs.end();
//		delIter--;
//		_recs.erase(delIter);
//		insert(dist, record, chromDir);
//		return true;
//	}
//
//}

//void RecDistList::insert(int dist, const Record *record, chromDirType chromDir) {
//	if ()
//	elemsType *elems = NULL;
//	distRecsType::iterator iter = _recs.find(dist);
//	if (iter == _recs.end()) {
//		//dist didn't exist before. Add it by
//		// creating new elems container, then
//		//putting this in the map.
//		elems = new elemsType();
//		_recs[dist] = elems;
//	} else {
//		elems = iter->second;
//	}
//	elems->push_back(pair<chromDirType, const Record *>(chromDir, record));
//	_totalRecs++;
//}

//if true, pos will be the idx the distance is at.
//if false, pos will be the idx to insert at.
bool RecDistList::find(int dist, int &pos) const {
	int lbound=0, ubound=_currNumIdxs-1, currVal =0;
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

int RecDistList::getMaxLeftEndPos() const {

	if (_empty) return -1;

	int maxDist =_distIndex[_currNumIdxs-1].first;
	const elemsType *elems = _allRecs[_distIndex[_currNumIdxs-1].second];
	for (int i=0; i < (int)elems->size(); i++) {
		if ((*elems)[i].first == LEFT) {
			return maxDist;
		}
	}
	return -1;

}


CloseSweep::CloseSweep(ContextClosest *context)
:	NewChromSweep(context),
 	_context(context),
 	_kClosest(_context->getNumClosestHitsWanted())
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
    return retVal;
 }

void CloseSweep::masterScan(RecordKeyVector &retList) {
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
    	const Record *cacheRec = cacheIter->value();
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


CloseSweep::rateOvlpType CloseSweep::considerRecord(const Record *cacheRec, int dbIdx, bool &stopScanning) {

	if (_context->diffNames() && cacheRec->getName() == _currQueryRec->getName()) {

		// We can ignore this, but we need to know whether to stop scanning.
		// do so IF:
		// 1 ) We are not ignoring downstream hits, AND
		// 2 ) The hit is after the query record, AND
		// 3 ) Some downstream hits have been found.
		if (!_context->ignoreDownstream() && cacheRec->after(_currQueryRec) && !_minDownstreamRecs[dbIdx]->empty()) {
			stopScanning = true;
		}

		// Secondly, we also want to know whether to delete this hit from the cache.
		//
		// TBD: Not sure how to determine this accurately. Leave it for now and hope the
		// performance doesn't suffer too badly.


		return IGNORE;
	}

	// If strand is specified, and either record has an unknown strand, ignore
	if ((_context->getSameStrand() || _context->getDiffStrand()) && ((_currQueryRec->getStrandVal() == Record::UNKNOWN) || cacheRec->getStrandVal() == Record::UNKNOWN)) {
		return IGNORE;
	}
	// If want same strand, and aren't sure they're the same, ignore
	if (_context->getSameStrand() &&  (_currQueryRec->getStrandVal() != cacheRec->getStrandVal())) {
		return IGNORE;
	}
	// If we want diff strand, and aren't sure they're different, ignore.
	if (_context->getDiffStrand() && (_currQueryRec->getStrandVal() == cacheRec->getStrandVal())) {
		return IGNORE;
	}

	// Now determine whether the hit and query intersect, and if so, what to do about it.
	int currDist = 0;

	if (intersects(_currQueryRec, cacheRec)) {

		// HIT INTERSECTS QUERY
		return tryToAddRecord(cacheRec, 0, dbIdx, stopScanning, OVERLAP, INTERSECT);

	} else if (cacheRec->after(_currQueryRec)) {

		// HIT IS TO THE RIGHT OF THE QUERY.

		 currDist = (cacheRec->getStartPos() - _currQueryRec->getEndPos()) + 1;
		 if (_context->signDistance()) {
			 if ((_context->getStrandedDistMode() == ContextClosest::A_DIST && _currQueryRec->getStrandVal() == Record::REVERSE) ||
				 (_context->getStrandedDistMode() == ContextClosest::B_DIST && cacheRec->getStrandVal() == Record::FORWARD))
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
			 if ((_context->getStrandedDistMode() == ContextClosest::A_DIST && _currQueryRec->getStrandVal() == Record::REVERSE) ||
				 (_context->getStrandedDistMode() == ContextClosest::B_DIST && cacheRec->getStrandVal() == Record::FORWARD))
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
	// This will control when the cahce is purged during the next query's sweep.
	setLeftClosestEndPos(dbIdx);



	RecDistList *upRecs = _minUpstreamRecs[dbIdx];
	RecDistList *downRecs = _minDownstreamRecs[dbIdx];
	RecDistList *overlaps = _overlapRecs[dbIdx];
	RecDistList::constIterType upIter = upRecs->begin();
	RecDistList::constIterType downIter = downRecs->begin();
	ContextClosest::tieModeType tieMode = _context->getTieMode();
	int upDist = INT_MAX;
	int downDist = INT_MAX;

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
		if (upDist < downDist || (tie && tieMode != ContextClosest::LAST_TIE)) {
			totalHitsUsed += addRecsToRetList(upRecs->allElems(upIter), 0 - upDist, retList);
			upIter++;
			usedUp = true;
		}
		if (downDist < upDist || (tie && tieMode != ContextClosest::FIRST_TIE)) {
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



int CloseSweep::addRecsToRetList(const RecDistList::elemsType *recs, int currDist, RecordKeyVector &retList) {

	int hitsUsed = 0;
	int numRecs = (int)recs->size(); //just to clean the code some.

	ContextClosest::tieModeType tieMode = _context->getTieMode();

	if (tieMode == ContextClosest::FIRST_TIE) {
		addSingleRec(recs->at(0).second, currDist, hitsUsed, retList);
		return 1;

	} else if (tieMode == ContextClosest::LAST_TIE) {
		addSingleRec(recs->at(numRecs-1).second, currDist, hitsUsed, retList);
		return 1;

	} else { //tieMode == ALL_TIES

		for (int i=0; i < numRecs; i++) {
			addSingleRec(recs->at(i).second, currDist, hitsUsed, retList);
		}
		return numRecs;
	}
}

void CloseSweep::addSingleRec(const Record *rec, int currDist, int &hitsUsed, RecordKeyVector &retList) {
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
	int i=0;
	for (RecordKeyVector::const_iterator_type iter = retList.begin(); iter != retList.end(); iter++) {
		int dist = _finalDistances[i];
		copyDists[i]._dist = abs(dist);
		copyDists[i]._rec = *iter;
		copyDists[i]._isNeg = dist < 0;
		i++;
	}
//	sort(copyDists.begin(), copyDists.end(), less<distanceTuple>());
	sort(copyDists.begin(), copyDists.end(), DistanceTupleSortAscFunctor());
	//now we want to build a map telling us what distances are tied,
	//and how many of each of these there are. Use a map<int, int>,
	//where the key is a distance (in absolute value) and the value
	//is the number of ties that that distance has.

	map<int, int> ties;
	for (i=0; i < numHits; i++) {
		int j = i +1;
		while ((j < numHits) && (copyDists[j]._dist == copyDists[i]._dist)) j++;
		if (j - i > 1) ties[copyDists[i]._dist] = j-i;
	}


	// Clear the original list and distances, and re-populate
	// until we have the desired number of hits, skipping
	// over any unwanted ties.
	retList.clearVector();
	_finalDistances.clear();

	int hitsUsed = 0;
	ContextClosest::tieModeType tieMode = _context->getTieMode();
	for (i=0; i < numHits && hitsUsed < _kClosest; i++) {
		int dist = copyDists[i]._dist;
		bool isNeg = copyDists[i]._isNeg;
		//see if this distance is tied with any other
		map<int, int>::iterator iter = ties.find(dist);
		if (iter != ties.end()) {
			//tie was found
			int numTies = iter->second;
			if (tieMode != ContextClosest::ALL_TIES) {
				if (tieMode == ContextClosest::FIRST_TIE) {
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
	const Record *dbRec = _currDbRecs[dbIdx];

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
		QuickString oldDbChrom(dbRec->getChrName());

		while (dbRec != NULL &&
				queryChromAfterDbRec(dbRec)) {
			_dbFRMs[dbIdx]->deleteRecord(dbRec);
			if (!nextRecord(false, dbIdx)) break;
			dbRec =  _currDbRecs[dbIdx];
			const QuickString &newDbChrom = dbRec->getChrName();
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
	RecDistList *upRecs = _minUpstreamRecs[dbIdx];
	RecDistList *downRecs = _minDownstreamRecs[dbIdx];

	int upDist = upRecs->getMaxLeftEndPos();
	int downDist = downRecs->getMaxLeftEndPos();

	if (upDist == -1 && downDist == -1) {
		//no records found to left of query,
		//can't set purge point, except for some special cases.
		if (purgePointException(dbIdx)) {
			_maxPrevLeftClosestEndPos[dbIdx] = max(_currQueryRec->getStartPos(), _maxPrevLeftClosestEndPos[dbIdx]);
		}
		return;
	}

	int leftMostEndPos = (_currQueryRec->getStartPos() - max(upDist, downDist));
	if ((!_context->getSameStrand() && !_context->getDiffStrand()) ||
		(_context->getSameStrand() && _currQueryRec->getStrandVal() == Record::FORWARD) ||
		(_context->getDiffStrand() && _currQueryRec->getStrandVal() == Record::REVERSE)) {

		_maxPrevLeftClosestEndPos[dbIdx] = max(leftMostEndPos, _maxPrevLeftClosestEndPos[dbIdx]);
	} else {
		_maxPrevLeftClosestEndPosReverse[dbIdx] = max(leftMostEndPos, _maxPrevLeftClosestEndPosReverse[dbIdx]);
	}
}

bool CloseSweep::beforeLeftClosestEndPos(int dbIdx, const Record *rec)
{
	int recEndPos = rec->getEndPos();
	int prevPos = _maxPrevLeftClosestEndPos[dbIdx];
	int prevPosReverse = _maxPrevLeftClosestEndPosReverse[dbIdx];

	if (!_context->getSameStrand() && !_context->getDiffStrand()) {
		return recEndPos < prevPos;
	} else {
		if (_context->getSameStrand()) {
			if (rec->getStrandVal() == Record::FORWARD) {
				return recEndPos < prevPos;
			} else {
				return recEndPos < prevPosReverse;
			}
		} else {
			//want diff strand
			if (rec->getStrandVal() == Record::FORWARD) {
				return recEndPos < prevPosReverse;
			} else {
				return recEndPos < prevPos;
			}
		}
	}
	return false;
}


void CloseSweep::clearClosestEndPos(int dbIdx)
{
	_maxPrevLeftClosestEndPos[dbIdx] = 0;
	_maxPrevLeftClosestEndPosReverse[dbIdx] = 0;
}

CloseSweep::rateOvlpType CloseSweep::tryToAddRecord(const Record *cacheRec, int dist, int dbIdx, bool &stopScanning, chromDirType chromDir, streamDirType streamDir) {

	//establish which set of hits we are looking at.
	RecDistList *useList = (streamDir == UPSTREAM ? _minUpstreamRecs[dbIdx] : (streamDir == INTERSECT ? _overlapRecs[dbIdx] : _minDownstreamRecs[dbIdx]));
	bool shouldIgnore = (streamDir == UPSTREAM ? _context->ignoreUpstream() : (streamDir == DOWNSTREAM ? _context->ignoreDownstream() :  _context->ignoreOverlaps()));

	if (chromDir == OVERLAP) {
		if (shouldIgnore || useList->addRec(0, cacheRec, RecDistList::OVERLAP)) {
		}
		return IGNORE;
	}

	if (chromDir == LEFT) {
		if (beforeLeftClosestEndPos(dbIdx, cacheRec)) {
			return DELETE;
		}
		if (!shouldIgnore) {
			useList->addRec(dist, cacheRec, RecDistList::LEFT);
		}
		return IGNORE;
	}
	// Here, chromDir == RIGHT
	else if (shouldIgnore || !useList->addRec(dist, cacheRec, RecDistList::RIGHT)) {
		// Stop scanning here, if we have no DIST mode.
		// If we do have a dist mode, or it's ref mode,
		// still stop scanning, UNLESS we only ignored the hit
		// because it and the query both have strands.
		if  (canStopScan(cacheRec, shouldIgnore, streamDir)) {
			stopScanning = true;
		} else {
//			printf("Couldn't stop scan.\n");
		}
	}
	return IGNORE;
}

bool CloseSweep::canStopScan(const Record *cacheRec, bool ignored, streamDirType streamDir) {
	if (!_context->hasStrandedDistMode() || (_context->getStrandedDistMode() == ContextClosest::REF_DIST) ||
			(cacheRec->getStrandVal() == Record::UNKNOWN || _currQueryRec->getStrandVal() == Record::UNKNOWN)) {
		return true;
	}
	//now we know we're in A_DIST or B_DIST mode, and the query and hit both have strands.
	if (!ignored){
		//ok, now we're only here because the record was out of range.
		return true;
	}
	// In A_DIST mode, when the query is negative, all records to the right are UPSTREAM.
	// pos query means all records are down stream. In these cases, it'll never stop scanning. Found that out the hard way.
	if ((_context->getStrandedDistMode() == ContextClosest::A_DIST) &&
		((streamDir == UPSTREAM && _currQueryRec->getStrandVal() == Record::REVERSE) ||
		 (streamDir == DOWNSTREAM && _currQueryRec->getStrandVal() == Record::FORWARD))) {
		return true;
	}

	return false;

}

bool CloseSweep::purgePointException(int dbIdx) {
	// we can only make an exception when ignoring upstream records, AND:
	// either have no dist mode, or it's REF dist mode, OR
	// either the query file or db file has unstranded records.
	return (_context->ignoreUpstream() &&
			((!_context->hasStrandedDistMode() || _context->getStrandedDistMode() == ContextClosest::REF_DIST) ||
			(!_context->getQueryFile()->recordsHaveStrand() || !_context->getDatabaseFile(dbIdx)->recordsHaveStrand())));
}
