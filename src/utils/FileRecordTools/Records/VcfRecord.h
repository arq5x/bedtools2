/*
 * VcfRecord.h
 *
 *  Created on: Apr 30, 2013
 *      Author: nek3d
 */

#ifndef VCFRECORD_H_
#define VCFRECORD_H_

#include "BedPlusInterval.h"

class SingleLineDelimTextFileReader;

class VcfRecord : public BedPlusInterval {
public:
	friend class FreeList<VcfRecord>;

	VcfRecord() {}
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::VCF_RECORD_TYPE; }

	//Note: using the assignment operator in a BedPlusInterval can potentially be a performance hit,
	//if the number of fields frequently differ between this object and the one being copied.
	const BedPlusInterval &operator=(const VcfRecord &other);

	virtual const QuickString &getField(int fieldNum) const;

protected:
	QuickString _varAlt;
	QuickString _varRef;
	void printOtherFields(QuickString &outBuf) const;
};

#endif /* VCFRECORD_H_ */
