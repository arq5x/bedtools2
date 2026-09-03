/*
 * BlockMgr.cpp
 *
 *  Created on: May 14, 2013
 *      Author: nek3d
 */

#include "BlockMgr.h"
#include "RecordMgr.h"
#include "Bed12Interval.h"
#include "BamRecord.h"
#include "ParseTools.h"
#include "api/BamAlignment.h"
#include "api/BamAux.h"

BlockMgr::BlockMgr(float overlapFractionA, float overlapFractionB, bool hasReciprocal)
: 	_blockRecordsMgr(NULL),
  	_breakOnDeletionOps(false),
  	_breakOnSkipOps(true),
  	_overlapFractionA(overlapFractionA),
  	_overlapFractionB(overlapFractionB),
  	_hasReciprocal(hasReciprocal)
{
	_blockRecordsMgr = new RecordMgr(_blockRecordsType);
}



BlockMgr::~BlockMgr()
{
	delete _blockRecordsMgr;
}

void BlockMgr::getBlocks(RecordKeyVector &keyList, bool &mustDelete)
{
	switch (keyList.getKey()->getType()) {
	case FileRecordTypeChecker::BED12_RECORD_TYPE:
		getBlocksFromBed12(keyList, mustDelete);
		break;

	case FileRecordTypeChecker::BAM_RECORD_TYPE:
		getBlocksFromBam(keyList, mustDelete);
		break;

	default:
		keyList.push_back(keyList.getKey());
		mustDelete = false;
		break;
	}
}

void BlockMgr::getBlocksFromBed12(RecordKeyVector &keyList, bool &mustDelete)
{
	const Bed12Interval *keyRecord = static_cast<const Bed12Interval *>(keyList.getKey());
	int blockCount = keyRecord->getBlockCount();

    if ( blockCount <= 0 ) {
    	mustDelete = false;
    	return;
    }

    int sizeCount = _blockSizeTokens.tokenize(keyRecord->getBlockSizes(), ',');
    int startCount = _blockStartTokens.tokenize(keyRecord->getBlockStarts(), ',');

    if (blockCount != sizeCount || sizeCount != startCount) {
    	fprintf(stderr, "Error: found wrong block counts while splitting entry.\n");
    	exit(-1);
    }

    for (int i=0; i < blockCount; i++) {
    	CHRPOS startPos = keyRecord->getStartPos() + str2chrPos(_blockStartTokens.getElem(i).c_str());
    	CHRPOS endPos = startPos + str2chrPos(_blockSizeTokens.getElem(i).c_str());

    	Record *record = allocateAndAssignRecord(keyRecord, startPos, endPos);
    	keyList.push_back(record);
    }
    mustDelete = true;
}

void BlockMgr::getBlocksFromBam(RecordKeyVector &keyList, bool &mustDelete)
{
	const BamRecord *keyRecord = static_cast<const BamRecord *>(keyList.getKey());
	const vector<BamTools::CigarOp> &cigarData = keyRecord->getCigarData();
	CHRPOS currPos = keyRecord->getStartPos();
	CHRPOS blockLength = 0;

	for (int i=0; i < (int)cigarData.size(); i++) {
		char opType = cigarData[i].Type;
		int opLen = (int)(cigarData[i].Length);

		switch(opType) {
		case 'I':
		case 'S':
		case 'P':
		case 'H':
			break;
		case 'M': case 'X': case '=':
			blockLength += opLen;
			break;
		case 'D':
		case 'N' :
			if ((opType == 'D' && !_breakOnDeletionOps) ||
					(opType == 'N' && !_breakOnSkipOps)) {
				blockLength += opLen;
			} else {
				keyList.push_back(allocateAndAssignRecord(keyRecord, currPos, currPos + blockLength));
				currPos += opLen + blockLength;
				blockLength = 0;
			}
			break;
		default:
			fprintf(stderr, "ERROR: Found invalid Cigar operation: %c.\n", opType);
			exit(1);
			break;
		}
	}
	if (blockLength > 0) {
		keyList.push_back(allocateAndAssignRecord(keyRecord, currPos, currPos + blockLength));
	}
	mustDelete = true;
}

Record *BlockMgr::allocateAndAssignRecord(const Record *keyRecord, CHRPOS startPos, CHRPOS endPos)
{
	Record *record = _blockRecordsMgr->allocateRecord();
	record->setChrName(keyRecord->getChrName());
	record->setChromId(keyRecord->getChromId());
	record->setStartPos(startPos);
	record->setEndPos(endPos);
	return record;
}

CHRPOS BlockMgr::getTotalBlockLength(RecordKeyVector &keyList) {
	CHRPOS sum = 0;
	for (RecordKeyVector::iterator_type iter = keyList.begin(); iter != keyList.end(); iter = keyList.next()) {
		const Record *record = *iter;
		sum += record->getEndPos() - record->getStartPos();
	}
	return sum;
}

void BlockMgr::deleteBlocks(RecordKeyVector &keyList)
{
	for (RecordKeyVector::iterator_type iter = keyList.begin(); iter != keyList.end(); iter = keyList.next()) {
		_blockRecordsMgr->deleteRecord(*iter);
	}
	keyList.clearVector();
}

int BlockMgr::getNonRedundantOverlap()
{
		// compute the non-redundant (merged) overlap
		// -f 0.2
	// ---------......--
    //         =      = 
    // Yes
	// 
	// ---------......--
    //                =
    //                =
    // No
	// ---------......--
    //         =......=
    // 
    // Yes
    //
	// ---------......--
    //                =
    //                 =
    // Yes
    //
	// ---------      
    //        ==
    //                
    // Yes
    //
	// ---------      
    //        =
    //        =
    //                
    // No	
    int totalNonRedundantOverlap = 0;
	if (_overlapStarts.size() == 1) 
	{
		CHRPOS start = _overlapStarts[0];
		CHRPOS end   = _overlapEnds[0];
		totalNonRedundantOverlap += end - start;
	}
	else if (_overlapStarts.size() > 1) {
		sort(_overlapStarts.begin(), _overlapStarts.end());
		sort(_overlapEnds.begin(), _overlapEnds.end());
		CHRPOS start = _overlapStarts[0];
		CHRPOS end   = _overlapEnds[0];
		for (int i = 1; i < _overlapStarts.size(); ++i)
		{
			CHRPOS currStart = _overlapStarts[i];
			CHRPOS currEnd = _overlapEnds[i];
			// new overlap block
			if (currStart >= end)
			{
				totalNonRedundantOverlap += end - start;
				start = currStart;
				end = currEnd;
			}
			// current overlap block continues
			else
			{
				start = currStart;
				end = currEnd;
			}
		}
		// account for last block
		totalNonRedundantOverlap += end - start;
	}
	return totalNonRedundantOverlap;
}

int BlockMgr::findBlockedOverlaps(RecordKeyVector &keyList, RecordKeyVector &hitList,
	                              RecordKeyVector &resultList, RecordKeyVector *overlapList)
{
	bool deleteKeyBlocks = false;
	if (keyList.empty()) {
		//get all the blocks for the query record, put them in it's list.
		getBlocks(keyList, deleteKeyBlocks);
	}
	_overlapBases.clear();
	_overlapStarts.clear();
	_overlapEnds.clear();
	CHRPOS keyBlocksSumLength = getTotalBlockLength(keyList);
	vector<CHRPOS> hitOverlapStarts;
	vector<CHRPOS> hitOverlapEnds;

	//Loop through every database record the query intersected with
	RecordKeyVector::iterator_type hitListIter = hitList.begin();
	for (; hitListIter != hitList.end(); hitListIter = hitList.next())
	{
		RecordKeyVector hitBlocks(*hitListIter);
		bool deleteHitBlocks = false;
		getBlocks(hitBlocks, deleteHitBlocks); //get all blocks for the hit record.
		CHRPOS hitBlockLength = getTotalBlockLength(hitBlocks); //total length of the blocks of this hit record.
		
		bool hitHasOverlap = false;
		CHRPOS totalDbOverlap = 0;
		hitOverlapStarts.clear();
		hitOverlapEnds.clear();
		//loop through every block of the database record.
		RecordKeyVector::iterator_type hitBlockIter = hitBlocks.begin();
		for (; hitBlockIter != hitBlocks.end(); hitBlockIter = hitBlocks.next()) 
		{
			//loop through every block of the query record.
			RecordKeyVector::iterator_type keyListIter = keyList.begin();
			for (; keyListIter != keyList.end(); keyListIter = keyList.next()) 
			{
				const Record *keyBlock = *keyListIter;
				const Record *hitBlock = *hitBlockIter;

				CHRPOS maxStart = max(keyBlock->getStartPos(), hitBlock->getStartPos());
				CHRPOS minEnd = min(keyBlock->getEndPos(), hitBlock->getEndPos());
				CHRPOS overlap  = minEnd - maxStart;
				if (overlap > 0) {
					hitHasOverlap = true;
					totalDbOverlap += overlap;
					hitOverlapStarts.push_back(maxStart);
					hitOverlapEnds.push_back(minEnd);
				}
			}
		}
		// -F and -r are tested per database record, as they are without -split:
		// the fraction of this record's blocks that the query's blocks cover.
		// (-f is cumulative over all records, below; see issue #750.)
		if (hitHasOverlap && (float)totalDbOverlap / (float)hitBlockLength < _overlapFractionB)
		{
			hitHasOverlap = false;
		}
		if (hitHasOverlap && _hasReciprocal &&
			(float)totalDbOverlap / (float)hitBlockLength < _overlapFractionA)
		{
			hitHasOverlap = false;
		}
		// add the database record to the list of potential (i.e., subject to -f)
		// hits if a sufficient overlap was observed
		if (hitHasOverlap)
		{
			resultList.push_back(*hitListIter);
			_overlapBases.push_back(totalDbOverlap);
			for (size_t i = 0; i < hitOverlapStarts.size(); i++) {
				_overlapStarts.push_back(hitOverlapStarts[i]);
				_overlapEnds.push_back(hitOverlapEnds[i]);
				if (overlapList != NULL) {
					overlapList->push_back(allocateAndAssignRecord(keyList.getKey(), hitOverlapStarts[i], hitOverlapEnds[i]));
				}
			}
		}
		if (deleteHitBlocks) {
			deleteBlocks(hitBlocks);
		}
	}

	int totalUniqueOverlap = getNonRedundantOverlap();

	// was there sufficient overlap with respect to -a? if not, delete.
	// This is the cumulative fraction of the query's blocks covered by
	// the records kept above (issue #750).
	if ((float) totalUniqueOverlap / (float)keyBlocksSumLength < _overlapFractionA) 
	{
		resultList.clearAll();
		_overlapBases.clear();
	}

	// clean up
	if (deleteKeyBlocks) {
		deleteBlocks(keyList);
	}
	resultList.setKey(keyList.getKey());
	return (int)resultList.size();
}

