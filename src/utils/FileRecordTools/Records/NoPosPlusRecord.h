/*
 * NoPosPlusRecord.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef NOPOSPLUSRECORD_H_
#define NOPOSPLUSRECORD_H_

#include "Record.h"
#include "PlusFields.h"

class SingleLineDelimTextFileReader;

class NoPosPlusRecord : public Record {
public:
	friend class FreeList<NoPosPlusRecord>;

	NoPosPlusRecord();
	virtual ~NoPosPlusRecord();
	bool initFromFile(FileReader *);
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void print(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return defaultNumFixedFields + (int)_plusFields.size(); }

	static bool isNumericField(int fieldNum);


protected:
	PlusFields _plusFields;
	static const int defaultNumFixedFields = 0;
};



#endif /* NOPOSPLUSRECORD_H_ */
