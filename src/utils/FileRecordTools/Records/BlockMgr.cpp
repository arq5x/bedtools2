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

BlockMgr::BlockMgr(float overlapFraction, bool hasReciprocal)
: 	_blockRecordsMgr(NULL),
  	_breakOnDeletionOps(false),
  	_breakOnSkipOps(true),
  	_overlapFraction(overlapFraction),
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
    	int startPos = keyRecord->getStartPos() + str2chrPos(_blockStartTokens.getElem(i).c_str());
    	int endPos = startPos + str2chrPos(_blockSizeTokens.getElem(i).c_str());

    	Record *record = allocateAndAssignRecord(keyRecord, startPos, endPos);
    	keyList.push_back(record);
    }
    mustDelete = true;
}

void BlockMgr::getBlocksFromBam(RecordKeyVector &keyList, bool &mustDelete)
{
	const BamRecord *keyRecord = static_cast<const BamRecord *>(keyList.getKey());
	const vector<BamTools::CigarOp> &cigarData = keyRecord->getCigarData();
	int currPos = keyRecord->getStartPos();
	int  blockLength = 0;

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

Record *BlockMgr::allocateAndAssignRecord(const Record *keyRecord, int startPos, int endPos)
{
	Record *record = _blockRecordsMgr->allocateRecord();
	record->setChrName(keyRecord->getChrName());
	record->setChromId(keyRecord->getChromId());
	record->setStartPos(startPos);
	record->setEndPos(endPos);
	return record;
}

int BlockMgr::getTotalBlockLength(RecordKeyVector &keyList) {
	int sum = 0;
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

int BlockMgr::findBlockedOverlaps(RecordKeyVector &keyList, RecordKeyVector &hitList,
	                              RecordKeyVector &resultList, RecordKeyVector *overlapList)
{
	bool deleteKeyBlocks = false;
	if (keyList.empty()) {
		//get all the blocks for the query record, put them in it's list.
		getBlocks(keyList, deleteKeyBlocks);
	}
	_overlapBases.clear();
	int keyBlocksSumLength = getTotalBlockLength(keyList);

	//Loop through every database record the query intersected with
	RecordKeyVector::iterator_type hitListIter = hitList.begin();
	for (; hitListIter != hitList.end(); hitListIter = hitList.next())
	{
		RecordKeyVector hitBlocks(*hitListIter);
		bool deleteHitBlocks = false;
		getBlocks(hitBlocks, deleteHitBlocks); //get all blocks for the hit record.
		int hitBlockSumLength = getTotalBlockLength(hitBlocks); //get total length of the bocks for the hitRecord.
		int totalHitOverlap = 0;
		bool hitHasOverlap = false;

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

				int maxStart = max(keyBlock->getStartPos(), hitBlock->getStartPos());
				int minEnd = min(keyBlock->getEndPos(), hitBlock->getEndPos());
				int overlap  = minEnd - maxStart;
				if (overlap > 0) {
					hitHasOverlap = true;
					if (overlapList != NULL) {
						overlapList->push_back(allocateAndAssignRecord(keyList.getKey(), maxStart, minEnd));
					}
					totalHitOverlap += overlap;
				}
			}
		}
		if (hitHasOverlap) {
			if ((float) totalHitOverlap / (float)keyBlocksSumLength >= _overlapFraction) {
				if (_hasReciprocal &&
						((float)totalHitOverlap / (float)hitBlockSumLength >= _overlapFraction)) {
					_overlapBases.push_back(totalHitOverlap);
					resultList.push_back(*hitListIter);
				} else if (!_hasReciprocal) {
					_overlapBases.push_back(totalHitOverlap);
					resultList.push_back(*hitListIter);
				}
			}
		}
		if (deleteHitBlocks) {
			deleteBlocks(hitBlocks);
		}
	}
	if (deleteKeyBlocks) {
		deleteBlocks(keyList);
	}
	resultList.setKey(keyList.getKey());
	return (int)resultList.size();
}

