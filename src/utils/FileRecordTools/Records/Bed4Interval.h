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
	virtual void print(string &outBuf) const;
	virtual void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	virtual void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED4_RECORD_TYPE; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 4; }
	static bool isNumericField(int fieldNum);


protected:
	virtual ~Bed4Interval();
};


#endif /* BED4INTERVAL_H_ */
