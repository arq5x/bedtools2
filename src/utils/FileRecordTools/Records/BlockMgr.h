/*
 * BlockMgr.h
 *
 *  Created on: May 14, 2013
 *      Author: nek3d
 */

#ifndef BLOCKMGR_H_
#define BLOCKMGR_H_

#include <vector>

//This class handles blocks inside of a larger record, such as BED12 and BAM records.
//Produce and manage seperate records for the sub-intervals inside the

#include "FileRecordTypeChecker.h"
#include "RecordKeyVector.h"

using namespace std;

class RecordMgr;

class BlockMgr {
public:
	BlockMgr(float overlapFractionA = 1E-9f, float overlapFractionB = 1E-9f, bool hasReciprocal = false);
	~BlockMgr();

	void getBlocks(RecordKeyVector &keyList, bool &mustDelete);
	
	int getNonRedundantOverlap();

	void deleteBlocks(RecordKeyVector &keyList);

	// Get the sum of the lengths of all blocks for a record.
	CHRPOS getTotalBlockLength(RecordKeyVector &keyList);


	// Determine which hits in the hitList intersect the hits in the keyList by comparing all blocks in each
	// and checking that their total intersection meets any overlapFraction and reciprocal criteria compared to
	// the total block lengths of the hitList and keyList. All hits that pass will be in the resultList.
	//
	// If useOverlappingSubBlocks is true, resultList will contain the sub-intervals of the hit blocks 
	// that overlap the query blocks.
	//
	// Return value is the number of hits in the result set.
	int findBlockedOverlaps(RecordKeyVector &keyList, RecordKeyVector &hitList,
							RecordKeyVector &resultList, RecordKeyVector *overlapList = NULL);

	// Same as above with no overlapList: the query's blocks come from key itself.
	int findBlockedOverlaps(const Record *key, RecordKeyVector &hitList, RecordKeyVector &resultList);

	//these are setting options for splitting BAM records
	void setBreakOnDeletionOps(bool val) { _breakOnDeletionOps = val; }
	void setBreakOnSkipOps(bool val) { _breakOnSkipOps = val; }
	int getOverlapBases(int hitIdx) const { return _overlapBases[hitIdx]; }

	Record *allocateAndAssignRecord(const Record *keyRecord, CHRPOS startPos, CHRPOS endPos);


private:
	RecordMgr *_blockRecordsMgr;
	bool _breakOnDeletionOps;
	bool _breakOnSkipOps;
	// keep track of the overlapping sub-intervals so
	// that we can tabulate non-redundant overlaps at the end
	vector<CHRPOS> _overlapBases;
	vector<CHRPOS> _overlapStarts;
	vector<CHRPOS> _overlapEnds;

	float _overlapFractionA;
	float _overlapFractionB;
	bool _hasReciprocal;
	vector<pair<size_t, size_t> > _blockSizeToks;
	vector<pair<size_t, size_t> > _blockStartToks;
	vector<pair<CHRPOS, CHRPOS> > _keyPairs;
	vector<pair<CHRPOS, CHRPOS> > _hitPairs;
	void getBlockPairs(const Record *rec, vector<pair<CHRPOS, CHRPOS> > &out);
	int finishBlockedOverlaps(const Record *key, CHRPOS keyBlocksSumLength,
			CHRPOS hitBlockSumLength, RecordKeyVector &resultList);

	// For now, all records will be split into Bed6 records.
	const static FileRecordTypeChecker::RECORD_TYPE _blockRecordsType = FileRecordTypeChecker::BED6_RECORD_TYPE;
	void getBlocksFromBed12(RecordKeyVector &keyList, bool &mustDelete);
	void getBlocksFromBam(RecordKeyVector &keyList, bool &mustDelete);

};


#endif /* BLOCKMGR_H_ */
