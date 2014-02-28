/*
 * Bed6Interval.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef BED6INTERVAL_H_
#define BED6INTERVAL_H_


#include "Bed3Interval.h"

class SingleLineDelimTextFileReader;

class Bed6Interval : public Bed3Interval {
public:
	friend class FreeList<Bed6Interval>;

	Bed6Interval();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED6_RECORD_TYPE; }

	virtual const QuickString &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 6; }
	static bool isNumericField(int fieldNum);


protected:
	virtual ~Bed6Interval();
};

#endif /* BED6INTERVAL_H_ */
