/*
 * GffRecord.h
 *
 *  Created on: Nov 13, 2012
 *      Author: nek3d
 */

#ifndef GFF_RECORD_H_
#define GFF_RECORD_H_

#include "Bed6Interval.h"

class SingleLineDelimTextFileReader;

class GffRecord : public Bed6Interval {
public:
	friend class FreeList<GffRecord>;

	GffRecord();
	virtual ~GffRecord();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::GFF_RECORD_TYPE; }
	virtual const QuickString &getSource() const { return _source; }
	virtual const QuickString &getFrame() const { return _frame; }
	virtual const QuickString &getGroup() const { return _group; }
	virtual int getNumFields() const { return _numFields; }
	virtual void setNumFields(int val) { _numFields = val; }

	virtual const QuickString &getField(int fieldNum) const;
	static bool isNumericField(int fieldNum);
	virtual bool isZeroBased() const {return false;};

protected:
	void printRemainingFields(QuickString &outbuf) const;

	int _numFields;
	QuickString _source;
	QuickString _frame;
	QuickString _group;

};



#endif /* GFF_RECORD_H_ */
