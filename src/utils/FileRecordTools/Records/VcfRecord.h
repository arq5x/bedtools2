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
	void print(QuickString &outBuf) const;
	void print(QuickString &outBuf, int start, int end) const;
	void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::VCF_RECORD_TYPE; }

	virtual const QuickString &getField(int fieldNum) const;
	static bool isNumericField(int fieldNum);

protected:
	QuickString _varRef;
	QuickString _varAlt;
	static const int numFixedFields = 6;

	void printOtherFields(QuickString &outBuf) const;
};

#endif /* VCFRECORD_H_ */
