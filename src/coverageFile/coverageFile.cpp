/*
 * coverageFile.cpp
 *
 *  Created on: May 8, 2015
 *      Author: nek3d
 */

#include "coverageFile.h"

CoverageFile::CoverageFile(ContextCoverage *context)
: IntersectFile(context),
 _depthArray(NULL),
 _depthArrayCapacity(0),
 _queryLen(0),
 _totalQueryLen(0),
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


void CoverageFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits) {
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
	IntersectFile::cleanupHits(hits);
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
		float depthPct = (float)basesAtDepth / (float)_totalQueryLen;

		_finalOutput = "all\t";
		_finalOutput.append(static_cast<uint32_t>(depth));
		_finalOutput.append("\t");
		_finalOutput.append(static_cast<uint32_t>(basesAtDepth));
		_finalOutput.append("\t");
		_finalOutput.append(static_cast<uint32_t>(_totalQueryLen));
		_finalOutput.append("\t");
		format(depthPct);

		outputMgr->printRecord(NULL, _finalOutput);
	}

}


void CoverageFile::makeDepthCount(RecordKeyVector &hits) {
	const Record *key = hits.getKey();
	_queryOffset = key->getStartPos();
	_queryLen = (size_t)(key->getEndPos() - _queryOffset);
	_totalQueryLen += _queryLen;

	//resize depth array if needed
	if (_depthArrayCapacity < _queryLen) {
		_depthArray = (size_t*)realloc(_depthArray, sizeof(size_t) * _queryLen);
		_depthArrayCapacity = _queryLen;
		memset(_depthArray, 0, sizeof(size_t) * _depthArrayCapacity);
	}

	//loop through hits, which may not be in sorted order, due to
	//potential multiple databases, and increment the depth array as needed.
	for (RecordKeyVector::const_iterator_type iter = hits.begin(); iter != hits.end(); iter = hits.next()) {
		const Record *dbRec = *iter;
		int dbStart = dbRec->getStartPos();
		int dbEnd = dbRec->getEndPos();
		int maxStart = max(_queryOffset, dbStart);
		int minEnd = min(dbEnd, key->getEndPos());

		for (int i=maxStart; i < minEnd; i++) {
			_depthArray[i - _queryOffset]++;
		}
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
	_finalOutput = static_cast<uint32_t>(hits.size());
	outputMgr->printRecord(hits.getKey(), _finalOutput);
}

void CoverageFile::doPerBase(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	//loop through all bases in query, printing full record and metrics for each
	const Record * queryRec = hits.getKey();
	for (size_t i= 0; i < _queryLen; i++) {
		_finalOutput = static_cast<uint32_t>(i+1);
		_finalOutput.append("\t");
		_finalOutput.append(static_cast<uint32_t>(_depthArray[i]));

		outputMgr->printRecord(queryRec, _finalOutput);
	}
}

void CoverageFile::doMean(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	size_t sum =0;
	for (size_t i= 0; i < _queryLen; i++) {
		sum += _depthArray[i];
	}
	format((float)sum / (float)_queryLen);
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
		float coveredBases = (float)numBasesAtDepth / (float)_queryLen;

		_finalOutput = static_cast<uint32_t>(depth);
		_finalOutput.append("\t");
		_finalOutput.append(static_cast<uint32_t>(numBasesAtDepth));
		_finalOutput.append("\t");
		_finalOutput.append(static_cast<uint32_t>(_queryLen));
		_finalOutput.append("\t");
		format(coveredBases);

		outputMgr->printRecord(hits.getKey(), _finalOutput);
	}

}

void CoverageFile::doDefault(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	size_t nonZeroBases = _queryLen - countBasesAtDepth(0);
	float coveredBases = (float)nonZeroBases / (float)_queryLen;

	_finalOutput = static_cast<uint32_t>(hits.size());
	_finalOutput.append("\t");
	_finalOutput.append(static_cast<uint32_t>(nonZeroBases));
	_finalOutput.append("\t");
	_finalOutput.append(static_cast<uint32_t>(_queryLen));
	_finalOutput.append("\t");
	format(coveredBases);

	outputMgr->printRecord(hits.getKey(), _finalOutput);
}

void CoverageFile::format(float val)
{
	memset(_floatValBuf, 0, floatValBufLen);
	sprintf(_floatValBuf, "%0.7f", val);
   _finalOutput.append(_floatValBuf);
}
