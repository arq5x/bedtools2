/*
 * FileRecordMergeMgr.cpp
 *
 *  Created on: Mar 19, 2014
 *      Author: nek3d
 */


#include "FileRecordMergeMgr.h"

FileRecordMergeMgr::FileRecordMergeMgr(const string & filename)
: FileRecordMgr(filename),
  _desiredStrand(ANY_STRAND),
  _maxDistance(0)
{
}

//Record *FileRecordMergeMgr::allocateAndGetNextMergedRecord(WANT_STRAND_TYPE desiredStrand, int maxDistance) {
//	RecordKeyVector recList;
//	if (!allocateAndGetNextMergedRecord(recList, desiredStrand, maxDistance)) {
//		return NULL;
//	}
//	deleteAllMergedItemsButKey(recList);
//	return const_cast<Record *>(recList.getKey()); //want key to be non-const
//}

Record *FileRecordMergeMgr::getNextRecord(RecordKeyVector *recList)
{
	//clear the recList if there is one, and if it has records
	// in it.
	if (recList != NULL && !recList->allClear()) {
		deleteMergedRecord(*recList);
	}

	_mustBeForward = _desiredStrand == SAME_STRAND_FORWARD;
	_mustBeReverse = _desiredStrand == SAME_STRAND_REVERSE;

	Record *startRecord = tryToTakeFromStorage();

	// if we couldn't use a previously stored record for starters,
	//then begin with a new one that matches strand criteria.
	while (startRecord == NULL) {
		startRecord = FileRecordMgr::getNextRecord();
		if (startRecord == NULL) { //hit EOF!!
			return NULL;
		}

		if ((_mustBeForward && (startRecord->getStrandVal() != Record::FORWARD)) || (_mustBeReverse && (startRecord->getStrandVal() != Record::REVERSE))) {
			//record is reverse, only want forward, OR record is forward, wanted reverse
			deleteRecord(startRecord);
			startRecord = NULL;
			continue;
		}
		if (startRecord->getStrandVal() == Record::UNKNOWN && _desiredStrand != ANY_STRAND) {
			//there is an unknown strand, but the user specified strandedness.
			deleteRecord(startRecord);
			startRecord = NULL;
		}
	}

	// OK!! We have a start record! Re-evaluate strand requirements for next recored.

	_mustBeForward = _desiredStrand == SAME_STRAND_FORWARD || (_desiredStrand == SAME_STRAND_EITHER && (startRecord->getStrandVal() == Record::FORWARD));
	_mustBeReverse = _desiredStrand == SAME_STRAND_REVERSE || (_desiredStrand == SAME_STRAND_EITHER && (startRecord->getStrandVal() == Record::REVERSE));
	bool mustKeepOpposite = (_desiredStrand == SAME_STRAND_EITHER);

	const string &currChrom = startRecord->getChrName();
	_foundChroms.insert(currChrom);

	bool madeComposite = false;
	if (recList != NULL) {
		recList->push_back(startRecord);
		recList->setKey(startRecord); //key of recList will just be the startRecord unless we're able to merge more.
	}

	Record::strandType currStrand = startRecord->getStrandVal();
	bool mustMatchStrand = _desiredStrand != ANY_STRAND;

	CHRPOS currEnd = startRecord->getEndPos();
	//now look for more records to merge with this one.
	//stop when they're out of range, not on the same chromosome, or we hit EOF.
	//ignore if they don't comply with strand.
	Record *nextRecord = NULL;
	while (nextRecord == NULL) {
		bool takenFromStorage = false;
		nextRecord = mustMatchStrand ? tryToTakeFromStorage(currStrand) : tryToTakeFromStorage();
		if (nextRecord == NULL) {
			nextRecord = FileRecordMgr::getNextRecord();
		} else {
			takenFromStorage = true;
		}
		if (nextRecord == NULL) { // EOF hit
			break;
		}
		//delete any record from file with an unknown strand if we are doing stranded merge, but first check
		//that it's chrom was the same and it's not out of range. If either is true, stop scanning.
		bool mustDelete = (mustMatchStrand && nextRecord->getStrandVal() == Record::UNKNOWN);

		//check that we are still on the same chromosome.
		const string &newChrom = nextRecord->getChrName();
		if (newChrom != currChrom) { //hit a different chromosome.
			//haven't seen this chromosome before, sort order is already enforced in the base class method.
			if (!mustDelete) {
				addToStorage(nextRecord);
			} else {
				deleteRecord(nextRecord);
			}
			nextRecord = NULL;
			break;
		}

		//check whether it's in range
		CHRPOS nextStart = nextRecord->getStartPos();
		if (nextStart > currEnd + _maxDistance) {
			//no, it's out of range.
			if (!mustDelete) {
				addToStorage(nextRecord);
			} else {
				deleteRecord(nextRecord);
			}
			nextRecord = NULL;
			break;
		}

		// NOW, going back, we can delete any unknown strand records. But don't stop scanning.
		if (mustDelete) {
			deleteRecord(nextRecord);
			nextRecord = NULL;
			continue;
		}
		//if taken from file, and wrong strand, store or delete.
		if (!takenFromStorage && ((_mustBeForward && (nextRecord->getStrandVal() != Record::FORWARD)) || (_mustBeReverse && (nextRecord->getStrandVal() != Record::REVERSE)))) {
			if (mustKeepOpposite) {
				addToStorage(nextRecord);
			} else {
				deleteRecord(nextRecord);
			}
			nextRecord = NULL;
			continue; //get the next record
		}
		//ok, they're on the same chrom and in range, and the strand is good. Do a merge.
		if (recList != NULL) recList->push_back(nextRecord);
		madeComposite = true;
		CHRPOS nextEnd = nextRecord->getEndPos();
		if (nextEnd > currEnd) {
			currEnd = nextEnd;
		}
		nextRecord = NULL;
	}
	if (madeComposite) {
		Record *newKey = _recordMgr->allocateRecord();
		(*newKey) = (*startRecord);
		newKey->setEndPos(currEnd);
		if (recList != NULL) recList->setKey(newKey);
		_totalMergedRecordLength += currEnd - newKey->getStartPos();
		return newKey;
	} else {
		_totalMergedRecordLength += currEnd - startRecord->getStartPos();
		return startRecord;
	}
//	_totalMergedRecordLength += (unsigned long)(recList->getKey()->getEndPos() - recList->getKey()->getStartPos());
//	return const_cast<Record *>(recList->getKey());
}

void FileRecordMergeMgr::addToStorage(Record *record) {
	//if the strand requirements are strict, and the record doesn't match,
	//store in the "round file".

	if ((_desiredStrand == SAME_STRAND_FORWARD && record->getStrandVal() != Record::FORWARD) ||
			(_desiredStrand == SAME_STRAND_REVERSE && record->getStrandVal() != Record::REVERSE) ||
			(_desiredStrand != ANY_STRAND && record->getStrandVal() == Record::UNKNOWN)) {
		deleteRecord(record);
		return;
	}
	_storedRecords.push(record);
}

Record *FileRecordMergeMgr::tryToTakeFromStorage() {
	Record *record = _storedRecords.top();
	if (record != NULL) {
		_storedRecords.pop();
	}
	return record;
}

Record *FileRecordMergeMgr::tryToTakeFromStorage(Record::strandType strand) {
	Record *record = _storedRecords.top(strand);
	if (record != NULL) {
		_storedRecords.pop(strand);
	}
	return record;
}

void FileRecordMergeMgr::deleteMergedRecord(RecordKeyVector &recList)
{
	deleteAllMergedItemsButKey(recList);
	deleteRecord(recList.getKey());
	recList.setKey(NULL);
}

bool FileRecordMergeMgr::eof(){
	return (_fileReader->eof() && _storedRecords.empty());
}


void FileRecordMergeMgr::deleteAllMergedItemsButKey(RecordKeyVector &recList) {
	//if the key is also in the list, this method won't delete it.
	for (RecordKeyVector::iterator_type iter = recList.begin(); iter != recList.end(); iter = recList.next()) {
		if (*iter == recList.getKey()) {
			continue;
		}
		deleteRecord(*iter);
	}
	recList.clearVector();
}


