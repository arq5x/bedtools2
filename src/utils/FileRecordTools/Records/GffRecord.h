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
	virtual void print(string &outBuf) const;
	virtual void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	virtual void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::GFF_RECORD_TYPE; }
	virtual const string &getSource() const { return _source; }
	virtual const string &getFrame() const { return _frame; }
	virtual const string &getGroup() const { return _group; }
	virtual int getNumFields() const { return _numFields; }
	virtual void setNumFields(int val) { _numFields = val; }

	virtual const string &getField(int fieldNum) const;
	static bool isNumericField(int fieldNum);
	virtual bool isZeroBased() const {return false;};

protected:
	void printRemainingFields(string &outbuf) const;

	int _numFields;
	string _source;
	string _frame;
	string _group;

};



#endif /* GFF_RECORD_H_ */
