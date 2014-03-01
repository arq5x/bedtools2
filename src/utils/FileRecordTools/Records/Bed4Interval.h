/*
 * Bed4Interval.h
 *
 *  Created on: Jun 27, 2013
 *      Author: nek3d
 */

#ifndef BED4INTERVAL_H_
#define BED4INTERVAL_H_



#include "Bed3Interval.h"

class SingleLineDelimTextFileReader;

class Bed4Interval : public Bed3Interval {
public:
	friend class FreeList<Bed4Interval>;

	Bed4Interval();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED4_RECORD_TYPE; }

	virtual const QuickString &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 4; }
	static bool isNumericField(int fieldNum);


protected:
	virtual ~Bed4Interval();
};


#endif /* BED4INTERVAL_H_ */
