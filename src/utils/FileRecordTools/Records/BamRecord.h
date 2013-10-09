/*
 * BamRecord.h
 *
 *  Created on: Dec 4, 2012
 *      Author: nek3d
 */

#ifndef BAMRECORD_H_
#define BAMRECORD_H_

#include "Bed6Interval.h"
#include "ParseTools.h"
#include "QuickString.h"
#include "api/BamAlignment.h"

class FileReader;
class BamFileReader;
class RecordKeyList;

class BamRecord : public Bed6Interval {
public:
	friend class FreeList<BamRecord>;

	BamRecord();
	virtual const BamRecord &operator=(const BamRecord &);
	bool initFromFile(FileReader *);
	virtual bool initFromFile(BamFileReader *);
	virtual void clear();
	virtual void print(QuickString &outBuf, int start, int end, RecordKeyList *keyList) const;
	virtual void print(QuickString &outBuf, RecordKeyList *keyList) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end, RecordKeyList *keyList) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual void printRemainingBamFields(QuickString &outBuf, RecordKeyList *keyList) const;
	virtual void printUnmapped(QuickString &outBuf) const;

	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BAM_RECORD_TYPE; }

	const BamTools::BamAlignment &getAlignment() const { return _bamAlignment; }
	int getBamChromId() const { return _bamChromId; }

protected:
	BamTools::BamAlignment _bamAlignment;
	int _bamChromId; //different from chromId, because BAM file may be in different order
	//than the genomeFile.

	virtual ~BamRecord();
	void printRemainingBamFields();

};


#endif /* BAMRECORD_H_ */
