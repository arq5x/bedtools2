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
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE; }

	virtual const QuickString &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 4; }

	static bool isNumericField(int fieldNum);

protected:
	virtual ~BedGraphInterval();
};

#endif /* BEDGRAPHINTERVAL_H_ */
