/*
 * RecordMgr.h
 *
 *  Created on: Jan 9, 2013
 *      Author: nek3d
 */

#ifndef RECORDMGR_H_
#define RECORDMGR_H_

//include headers for all Records and derivative classes
#include "Record.h"
#include "EmptyRecord.h"
#include "Bed3Interval.h"
#include "Bed4Interval.h"
#include "Bed5Interval.h"
#include "BedGraphInterval.h"
#include "Bed6Interval.h"
#include "BedPlusInterval.h"
#include "Bed12Interval.h"
#include "BamRecord.h"
#include "VcfRecord.h"
#include "GffRecord.h"
#include "GffPlusRecord.h"
#include "NoPosPlusRecord.h"

#include "FileRecordTypeChecker.h"

class Record;

class RecordMgr {
public:
	RecordMgr(FileRecordTypeChecker::RECORD_TYPE recType, int blockSize = 512);
	~RecordMgr();

	Record *allocateRecord();

	//WARNING: Even though the deleteRecord method takes a const pointer, user should only pass objects
	//that are safe to delete! Const-ness will be cast away internally.
	void deleteRecord(const Record *);

private:
	FileRecordTypeChecker::RECORD_TYPE _recordType;
	void *_freeList;
	int _freeListBlockSize;
};



#endif /* RECORDMGR_H_ */
