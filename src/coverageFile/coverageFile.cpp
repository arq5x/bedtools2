/*
 * coverageFile.cpp
 *
 *  Created on: May 8, 2015
 *      Author: nek3d
 */

#include "coverageFile.h"
#include <iomanip>

CoverageFile::CoverageFile(ContextCoverage *context)
: IntersectFile(context),
 _depthArray(NULL),
 _depthArrayCapacity(0),
 _queryLen(0),
 _totalQueryLen(0),
 _hitCount(0),
 _queryOffset(0),
 _floatValBuf(NULL)
{
	//allocate and initialize depth array.
	//Use C's malloc and free becauase we want
	//to be able to realloc.
	//Not using vector because arrays are faster.
	size_t newSize = sizeof(size_t) * DEFAULT_DEPTH_CAPACITY;
	_depthArray = (size_t *)malloc(newSize);
	memset(_depthArray, 0, newSize);
	_depthArrayCapacity = DEFAULT_DEPTH_CAPACITY;

	_floatValBuf = new char[floatValBufLen];
}

CoverageFile::~CoverageFile() {
	free(_depthArray);
	delete _floatValBuf;
}


void CoverageFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits) 
{
	makeDepthCount(hits);
	_finalOutput.clear();

	switch(upCast(_context)->getCoverageType()) {
	case ContextCoverage::COUNT:
	   doCounts(outputMgr, hits);
	   break;

	case ContextCoverage::PER_BASE:
	   doPerBase(outputMgr, hits);
	   break;

	case ContextCoverage::MEAN:
	   doMean(outputMgr, hits);
	   break;

	case ContextCoverage::HIST:
	   doHist(outputMgr, hits);
	   break;

	case ContextCoverage::DEFAULT:
	default:
	   doDefault(outputMgr, hits);
	   break;
	}
}

void CoverageFile::cleanupHits(RecordKeyVector &hits) {
	if (upCast(_context)->getObeySplits()) {
		upCast(_context)->getSplitBlockInfo()->deleteBlocks(hits);
	} else {
		IntersectFile::cleanupHits(hits);
	}
	memset(_depthArray, 0, sizeof(size_t) * _queryLen);
}

void CoverageFile::giveFinalReport(RecordOutputMgr *outputMgr) {

	//only give report for histogram option
	if (upCast(_context)->getCoverageType() != ContextCoverage::HIST) {
		return;
	}

	
	for (depthMapType::iterator iter = _finalDepthMap.begin(); iter != _finalDepthMap.end(); iter++) {
		size_t depth = iter->first;
		size_t basesAtDepth = iter->second;
		//cout << "x\n";
		float depthPct = (float)basesAtDepth / (float)_totalQueryLen;
		//cout << "y\n";
		ostringstream s;
		s << "all\t";
		s << depth;
		s << "\t";
		s << basesAtDepth;
		s << "\t";
		s << _totalQueryLen;
		s << "\t";
		char *depthPctString;
		asprintf(&depthPctString, "%0.7f", depthPct);
		s << depthPctString;
		_finalOutput = s.str();

		outputMgr->printRecord(NULL, _finalOutput);
	}
}

void CoverageFile::makeDepthCount(RecordKeyVector &hits) {
	const Record *key = hits.getKey();
	_queryOffset = key->getStartPos();
	_queryLen = (size_t)(key->getEndPos() - _queryOffset);
	_totalQueryLen += _queryLen;

	// resize depth array if needed
	if (_depthArrayCapacity < _queryLen) {
		_depthArray = (size_t*)realloc(_depthArray, sizeof(size_t) * _queryLen);
		_depthArrayCapacity = _queryLen;
		memset(_depthArray, 0, sizeof(size_t) * _depthArrayCapacity);
	}
	_hitCount = 0;
	// Here we do not want extra logic for splits, b/c they're already handled in how we're getting hits.
	//loop through hits, which may not be in sorted order, due to
	//potential multiple databases, and increment the depth array as needed.
	for (RecordKeyVector::iterator_type iter = hits.begin(); iter != hits.end(); iter = hits.next())
	{
		const Record *dbRec = *iter;
		CHRPOS dbStart = dbRec->getStartPos();
		CHRPOS dbEnd = dbRec->getEndPos();
		CHRPOS maxStart = max(_queryOffset, dbStart);
		CHRPOS minEnd = min(dbEnd, key->getEndPos());

		for (CHRPOS i=maxStart; i < minEnd; i++) {
			_depthArray[i - _queryOffset]++;
		}
		_hitCount++;
	}
}

size_t CoverageFile::countBasesAtDepth(size_t depth) {
	size_t retCount = 0;
	for (size_t i=0; i < _queryLen; i++) {
		if (_depthArray[i] == depth) {
			retCount++;
		}
	}
	return retCount;
}

void CoverageFile::doCounts(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	ostringstream s;
	s << _hitCount;
	_finalOutput.append(s.str());
	outputMgr->printRecord(hits.getKey(), _finalOutput);
}

void CoverageFile::doPerBase(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	//loop through all bases in query, printing full record and metrics for each
	
	Record * queryRec = hits.getKey();
	for (size_t i= 0; i < _queryLen; i++) {
		ostringstream s;
		s << (i+1);
		s << "\t";
		s << _depthArray[i];
		_finalOutput = s.str();
		outputMgr->printRecord(queryRec, _finalOutput);
	}
}

void CoverageFile::doMean(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	size_t sum =0;
	for (size_t i= 0; i < _queryLen; i++) {
		sum += _depthArray[i];
	}
	ostringstream s;
	float mean = ((float)sum / (float)_queryLen);
	char *meanString;
	asprintf(&meanString, "%0.7f", mean);
	s << meanString;
	_finalOutput.append(s.str());
	outputMgr->printRecord(hits.getKey(), _finalOutput);
}


void CoverageFile::doHist(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	//make a map of depths to num bases with that depth
	_currDepthMap.clear();
	for (size_t i=0; i < _queryLen; i++) {
		_currDepthMap[_depthArray[i]]++;
		_finalDepthMap[_depthArray[i]]++;
	}

	for (depthMapType::iterator iter = _currDepthMap.begin(); iter != _currDepthMap.end(); iter++) {
		size_t depth = iter->first;
		size_t numBasesAtDepth = iter->second;
		float coveredFraction = (float)numBasesAtDepth / (float)_queryLen;

		ostringstream s;
		s << depth;
		s << "\t";
		s << numBasesAtDepth;
		s << "\t";
		s << _queryLen;
		s << "\t";
		char *coveredFractionString;
		asprintf(&coveredFractionString, "%0.7f", coveredFraction);
		s << coveredFractionString;
		_finalOutput = s.str();
		outputMgr->printRecord(hits.getKey(), _finalOutput);
	}
}

void CoverageFile::doDefault(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	size_t nonZeroBases = _queryLen - countBasesAtDepth(0);
	float coveredFraction = (float)nonZeroBases / (float)_queryLen;

	ostringstream s;
	s << _hitCount;
	s << "\t";
	s << nonZeroBases;
	s << "\t";
	s << _queryLen;
	s << "\t";
	char *coveredFractionString;
	asprintf(&coveredFractionString, "%0.7f", coveredFraction);
	s << coveredFractionString;
	_finalOutput = s.str();
	outputMgr->printRecord(hits.getKey(), _finalOutput);
}

void CoverageFile::checkSplits(RecordKeyVector &hitSet)
{
	// When using coverage, we need a list of the sub-intervals of coverage
	// so that per-base depth can be properly calculated when obeying splits
	if (upCast(_context)->getObeySplits()) {
		RecordKeyVector keySet(hitSet.getKey());
		RecordKeyVector resultSet(hitSet.getKey());
		RecordKeyVector overlapSet(hitSet.getKey());
		upCast(_context)->getSplitBlockInfo()->findBlockedOverlaps(keySet, hitSet, resultSet, &overlapSet);
		//hitSet.clearAll();
		hitSet.swap(overlapSet);
	}
}
