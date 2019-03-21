/*
 * Bed5Interval.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef BED5INTERVAL_H_
#define BED5INTERVAL_H_


#include "Bed3Interval.h"

class SingleLineDelimTextFileReader;

class Bed5Interval : public Bed3Interval {
public:
	friend class FreeList<Bed5Interval>;

	Bed5Interval();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void print(string &outBuf) const;
	virtual void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	virtual void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED5_RECORD_TYPE; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 5; }
	static bool isNumericField(int fieldNum);


protected:
	virtual ~Bed5Interval();
};

#endif /* BED5INTERVAL_H_ */
