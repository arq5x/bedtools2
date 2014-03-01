/*
 * BedPlusInterval.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef BEDPLUSINTERVAL_H_
#define BEDPLUSINTERVAL_H_

#include "Bed6Interval.h"
#include <vector>

class SingleLineDelimTextFileReader;

class BedPlusInterval : public Bed6Interval {
public:
	friend class FreeList<BedPlusInterval>;

	BedPlusInterval();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED_PLUS_RECORD_TYPE; }

	//Note: using the assignment operator in a BedPlusInterval can potentially be a performance hit,
	//if the number of fields frequently differ between this object and the one being copied.
	const BedPlusInterval &operator=(const BedPlusInterval &other);

	virtual const QuickString &getField(int fieldNum) const;
	virtual int getNumFields() const  { return startOtherIdx + _otherIdxs.size(); }

	virtual void setField(int fieldNum, const QuickString &str) { (*(_otherIdxs[fieldNum])) = str; }
	virtual void setField(int fieldNum, const string &str) { (*(_otherIdxs[fieldNum])) = str; }
	virtual void setField(int fieldNum, const char *str) { (*(_otherIdxs[fieldNum])) = str; }
	virtual void setNumPrintFields(int num) { _numPrintFields = num; }
	virtual int getNumPrintFields() const { return _numPrintFields; }
	static bool isNumericField(int fieldNum);


protected:
	virtual ~BedPlusInterval();
	bool initOtherFieldsFromFile(SingleLineDelimTextFileReader *fileReader);


	vector<QuickString *> _otherIdxs;
	static const int startOtherIdx = 6; //first six fields have names, and are not stored in otherIdxs.
	int _numPrintFields;
};



#endif /* BEDPLUSINTERVAL_H_ */
