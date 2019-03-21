/*
 * GffRecord.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef GFFPLUSRECORD_H_
#define GFFPLUSRECORD_H_

#include "GffRecord.h"
#include "PlusFields.h"

class SingleLineDelimTextFileReader;

class GffPlusRecord : public GffRecord {
public:
	friend class FreeList<GffPlusRecord>;

	GffPlusRecord();
	virtual ~GffPlusRecord();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void print(string &outBuf) const;
	virtual void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	virtual void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::GFF_PLUS_RECORD_TYPE; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return (int)GffRecord::getNumFields() + (int)_plusFields.size(); }

	virtual void setNumPrintFields(int num) { _numPrintFields = num; }
	virtual int getNumPrintFields() const { return _numPrintFields; }
	static bool isNumericField(int fieldNum);
	virtual bool isZeroBased() const {return false;};

protected:
	PlusFields _plusFields;
	int _numPrintFields;

};



#endif /* GFFPLUSRECORD_H_ */
