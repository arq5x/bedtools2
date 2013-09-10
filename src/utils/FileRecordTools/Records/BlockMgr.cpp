/*
 * BlockMgr.cpp
 *
 *  Created on: May 14, 2013
 *      Author: nek3d
 */

#include "BlockMgr.h"
#include "RecordMgr.h"
#include "Context.h"
#include "Bed12Interval.h"
#include "BamRecord.h"
#include "ParseTools.h"
#include "api/BamAlignment.h"
#include "api/BamAux.h"

BlockMgr::BlockMgr()
: 	_blockRecordsMgr(NULL),
  	_breakOnDeletionOps(false),
  	_breakOnSkipOps(true)
{
	_blockRecordsMgr = new RecordMgr(_blockRecordsType);
}



BlockMgr::~BlockMgr()
{
	delete _blockRecordsMgr;
}

void BlockMgr::getBlocks(RecordKeyList &keyList, bool &mustDelete)
{
	switch (keyList.getKey()->getType()) {
	case FileRecordTypeChecker::BED12_RECORD_TYPE:
		getBlocksFromBed12(keyList, mustDelete);
		break; //not necessary after return, but here in case code is later modified.

	case FileRecordTypeChecker::BAM_RECORD_TYPE:
		getBlocksFromBam(keyList, mustDelete);
		break;

	default:
		keyList.push_back(keyList.getKey());
		mustDelete = false;
		break;
	}
}



void BlockMgr::getBlocksFromBed12(RecordKeyList &keyList, bool &mustDelete)
{
	const Bed12Interval *keyRecord = static_cast<const Bed12Interval *>(keyList.getKey());
	int blockCount = keyRecord->getBlockCount();

    if ( blockCount <= 0 ) {
    	mustDelete = false;
    	return;
    }

    vector<QuickString> sizes;
    vector<QuickString> starts;

    int sizeCount = Tokenize(keyRecord->getBlockSizes(), sizes, ',', blockCount);
    int startCount = Tokenize(keyRecord->getBlockStarts(), starts, ',', blockCount);

    if (blockCount != sizeCount || sizeCount != startCount) {
    	fprintf(stderr, "Error: found wrong block counts while splitting entry.\n");
    	exit(-1);
    }

    for (int i=0; i < blockCount; i++) {
    	int startPos = keyRecord->getStartPos() + str2chrPos(starts[i].c_str());
    	int endPos = startPos + str2chrPos(sizes[i].c_str());

    	const Record *record = allocateAndAssignRecord(keyRecord, startPos, endPos);
    	keyList.push_back(record);
    }
    mustDelete = true;
}

void BlockMgr::getBlocksFromBam(RecordKeyList &keyList, bool &mustDelete)
{
	const BamRecord *keyRecord = static_cast<const BamRecord *>(keyList.getKey());
	const vector<BamTools::CigarOp> &cigarData = keyRecord->getAlignment().CigarData;

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
		case 'M':
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
	record->setName(keyRecord->getName());
	record->setScore(keyRecord->getScore());
	record->setStrand(keyRecord->getStrand());

	return record;
}

int BlockMgr::getTotalBlockLength(RecordKeyList &keyList) {
	int sum = 0;
	for (RecordKeyList::const_iterator_type iter = keyList.begin(); iter != keyList.end(); iter = keyList.next()) {
		const Record *record = iter->value();
		sum += record->getEndPos() - record->getStartPos();
	}
	return sum;
}

void BlockMgr::deleteBlocks(RecordKeyList &keyList)
{
	for (RecordKeyList::const_iterator_type iter = keyList.begin(); iter != keyList.end(); iter = keyList.next()) {
		_blockRecordsMgr->deleteRecord(iter->value());
	}
	keyList.clearList();

}


int BlockMgr::findBlockedOverlaps(RecordKeyList &keyList, RecordKeyList &hitList, RecordKeyList &resultList)
{
	bool deleteKeyBlocks = false;
	if (keyList.empty()) {
		getBlocks(keyList, deleteKeyBlocks);
	}
	int keyBlocksSumLength = getTotalBlockLength(keyList);
	for (RecordKeyList::const_iterator_type hitListIter = hitList.begin(); hitListIter != hitList.end(); hitListIter = hitList.next()) {
		RecordKeyList hitBlocks(hitListIter->value());
		bool deleteHitBlocks = false;
		getBlocks(hitBlocks, deleteHitBlocks);
		int hitBlockSumLength = getTotalBlockLength(hitBlocks);
		int totalHitOverlap = 0;
		bool hitHasOverlap = false;
		for (RecordKeyList::const_iterator_type hitBlockIter = hitBlocks.begin(); hitBlockIter != hitBlocks.end(); hitBlockIter = hitBlocks.next()) {
			for (RecordKeyList::const_iterator_type keyListIter = keyList.begin(); keyListIter != keyList.end(); keyListIter = keyList.next()) {
				const Record *keyBlock = keyListIter->value();
				const Record *hitBlock = hitBlockIter->value();

				int maxStart = max(keyBlock->getStartPos(), hitBlock->getStartPos());
				int minEnd = min(keyBlock->getEndPos(), hitBlock->getEndPos());
				int overlap  = minEnd - maxStart;
				if (overlap > 0) {
					hitHasOverlap = true;
					totalHitOverlap += overlap;
				}

			}
		}
		if (hitHasOverlap) {
			if ((float) totalHitOverlap / (float)keyBlocksSumLength > _context->getOverlapFraction()) {
				if (_context->getReciprocal() &&
						((float)totalHitOverlap / (float)hitBlockSumLength > _context->getOverlapFraction())) {
					resultList.push_back(hitListIter->value());
				} else if (!_context->getReciprocal()) {
					resultList.push_back(hitListIter->value());
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

