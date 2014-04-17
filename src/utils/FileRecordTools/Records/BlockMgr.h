/*
 * BlockMgr.h
 *
 *  Created on: May 14, 2013
 *      Author: nek3d
 */

#ifndef BLOCKMGR_H_
#define BLOCKMGR_H_

//This class handles blocks inside of a larger record, such as BED12 and BAM records.
//Produce and manage seperate records for the sub-intervals inside the

using namespace std;

#include "FileRecordTypeChecker.h"
#include "RecordKeyList.h"


class RecordMgr;

class BlockMgr {
public:
	BlockMgr(float overlapFraction = 1E-9, bool hasReciprocal = false);
	~BlockMgr();

	// Return value is the number of blocks this main record has been split into.
	void getBlocks(RecordKeyList &keyList, bool &mustDelete);
	void deleteBlocks(RecordKeyList &keyList);

	// Get the sum of the lengths of all blocks for a record.
	int getTotalBlockLength(RecordKeyList &keyList);

	// Determine which hits in the hitList intersect the hits in the keyList by comparing all blocks in each
	// and checking that their total intersection meets any overlapFraction and reciprocal criteria compared to
	// the total block lengths of the hitList and keyList. All hits that pass will be in the resultList.
	// Return value is the number of hits in the result set.

	int findBlockedOverlaps(RecordKeyList &keyList, RecordKeyList &hitList, RecordKeyList &resultList);

	//these are setting options for splitting BAM records
	void setBreakOnDeletionOps(bool val) { _breakOnDeletionOps = val; }
	void setBreakOnSkipOps(bool val) { _breakOnSkipOps = val; }
	int getOverlapBases(int hitIdx) const { return _overlapBases[hitIdx]; }

private:
	RecordMgr *_blockRecordsMgr;
	bool _breakOnDeletionOps;
	bool _breakOnSkipOps;
	vector<int> _overlapBases;

	float _overlapFraction;
	bool _hasReciprocal;
	Tokenizer _blockSizeTokens;
	Tokenizer _blockStartTokens;

	// For now, all records will be split into Bed6 records.
	const static FileRecordTypeChecker::RECORD_TYPE _blockRecordsType = FileRecordTypeChecker::BED6_RECORD_TYPE;
	void getBlocksFromBed12(RecordKeyList &keyList, bool &mustDelete);
	void getBlocksFromBam(RecordKeyList &keyList, bool &mustDelete);

	Record *allocateAndAssignRecord(const Record *keyRecord, int startPos, int endPos);


};


#endif /* BLOCKMGR_H_ */
