/*
 * Bed6Interval.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef BEDGRAPHINTERVAL_H_
#define BEDGRAPHINTERVAL_H_


#include "Bed3Interval.h"

class SingleLineDelimTextFileReader;

class BedGraphInterval : public Bed3Interval {
public:
	friend class FreeList<BedGraphInterval>;

	BedGraphInterval();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void print(string &outBuf) const;
	virtual void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	virtual void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 4; }

	static bool isNumericField(int fieldNum);

protected:
	virtual ~BedGraphInterval();
};

#endif /* BEDGRAPHINTERVAL_H_ */
