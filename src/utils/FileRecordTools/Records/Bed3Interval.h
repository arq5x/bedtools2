/*
 * Bed3Interval.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef BED3INTERVAL_H_
#define BED3INTERVAL_H_

#include "Record.h"
#include "ParseTools.h"
#include "string.h"

class FileReader;
class SingleLineDelimTextFileReader;


class Bed3Interval : public Record {
public:
	friend class FreeList<Bed3Interval>;

	Bed3Interval();
	bool initFromFile(FileReader *);
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void print(string &outBuf) const;
	virtual void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	virtual void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED3_RECORD_TYPE; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 3; }

	static bool isNumericField(int fieldNum);
	virtual ~Bed3Interval();

protected:

	bool _skipFirstThreeFieldsWhenPrinting;

};


#endif /* BED3INTERVAL_H_ */
