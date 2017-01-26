/*
 * FileRecordMergeMgr.h
 *
 *  Created on: Mar 19, 2014
 *      Author: nek3d
 */

#ifndef FILERECORDMERGEMGR_H_
#define FILERECORDMERGEMGR_H_

#include "FileRecordMgr.h"
#include "StrandQueue.h"

class FileRecordMergeMgr : public FileRecordMgr {

public:
	FileRecordMergeMgr(const string & filename);

	//////////////////////////////////////////////////////////////////////////////////
	//
	// 			MERGED RECORDS
	//
	// This will give a single "meta" record containing "flattened" or merged records.
	//
	// Pass an empty RecordKeyVector. When done, will have a pair: 1st is the final merged record,
	//			second is list of constituent Records merged.
	//
	///////////////////////////////////////////////////////////////////////////////////

	Record *getNextRecord(RecordKeyVector *keyList = NULL);
	void deleteMergedRecord(RecordKeyVector &recList); // MUST use this method for cleanup!

	bool eof();


	typedef enum { SAME_STRAND_FORWARD, //must all be forward strand
			SAME_STRAND_REVERSE, //must all be reverse strand
			SAME_STRAND_EITHER, //must be same strand, but can be either forward or reverse
			ANY_STRAND } //do no care about strand (Default value)
	WANTED_STRAND_TYPE;

	void setStrandType(WANTED_STRAND_TYPE strand) { _desiredStrand = strand; }
	void setMaxDistance(int maxDistance) { _maxDistance = maxDistance; }

private:

	WANTED_STRAND_TYPE _desiredStrand;
	int _maxDistance;
	StrandQueue _storedRecords;

	void deleteAllMergedItemsButKey(RecordKeyVector &recList);
	void addToStorage(Record *record);
	Record *tryToTakeFromStorage();
	Record *tryToTakeFromStorage(Record::strandType strand);
};


#endif /* FILERECORDMERGEMGR_H_ */
