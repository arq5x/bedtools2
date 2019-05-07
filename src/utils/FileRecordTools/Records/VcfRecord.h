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
	void print(string &outBuf) const;
	void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::VCF_RECORD_TYPE; }

	virtual const string &getField(int fieldNum) const;
	static bool isNumericField(int fieldNum);

	virtual bool isZeroBased() const {return false;};
protected:
	string _varRef;
	string _varAlt;
	static const int numFixedFields = 6;

	void printOtherFields(string &outBuf) const;
};

#endif /* VCFRECORD_H_ */
