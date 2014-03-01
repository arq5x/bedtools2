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
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED5_RECORD_TYPE; }

	virtual const QuickString &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 5; }
	static bool isNumericField(int fieldNum);


protected:
	virtual ~Bed5Interval();
};

#endif /* BED5INTERVAL_H_ */
