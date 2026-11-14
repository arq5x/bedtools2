/*
 * Bed12Interval.h
 *
 *  Created on: Nov 28, 2012
 *      Author: nek3d
 */

#ifndef BED12INTERVAL_H_
#define BED12INTERVAL_H_

#include "Bed6Interval.h"
#include "lineFileUtilities.h"
#include <numeric>
class SingleLineDelimTextFileReader;


class Bed12Interval : public Bed6Interval {
public:
	friend class FreeList<Bed12Interval>;

	Bed12Interval();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void print(string &outBuf) const;
	virtual void print(string &outBuf, CHRPOS start, CHRPOS end) const;
	virtual void print(string &outBuf, const string & start, const string & end) const;
	virtual void printNull(string &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED12_RECORD_TYPE; }

	const Bed12Interval &operator=(const Bed12Interval &other);

	virtual CHRPOS getThickStart() const { return _thickStart; }
	virtual void setThickStart(CHRPOS thickStart)  { _thickStart = thickStart; _lineValid = false; }

	virtual CHRPOS getThickEnd() const { return _thickEnd; }
	virtual void setThickEnd(CHRPOS thickEnd) { _thickEnd = thickEnd; _lineValid = false; }

	virtual const string & getItemRGB() const { return _itemRGB; }
	virtual void setItemRGB(const string & rgb) { _itemRGB = rgb; _lineValid = false; }
	virtual void setItemRGB(const char *rgb) { _itemRGB = rgb; _lineValid = false; }

	virtual int getBlockCount() const { return _blockCount; }
	virtual void setBlockCount(int blockCount) { _blockCount = blockCount; _lineValid = false; }

	virtual const string & getBlockSizes() const { return _blockSizes; }
	virtual void setBlockSizes(const string & blockSizes) { _blockSizes = blockSizes; _lineValid = false; }
	virtual void setBlockSizes(const char *blockSizes) { _blockSizes = blockSizes; _lineValid = false; }

	virtual const string & getBlockStarts() const { return _blockStarts; }
	virtual void setBlockStarts(const string & blockStarts) { _blockStarts = blockStarts; _lineValid = false; }
	virtual void setBlockStarts(const char *blockStarts) { _blockStarts = blockStarts; _lineValid = false; }

	virtual const string &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 12; }
	static bool isNumericField(int fieldNum);
	CHRPOS getLength(bool obeySplits) const;


protected:
	virtual ~Bed12Interval();

	CHRPOS _thickStart;
	CHRPOS _thickEnd;
	string _itemRGB;
	int _blockCount;
	string _blockSizes;
	string _blockStarts;

	//store int values as their originial strings for faster outputting.
	string _thickStartStr;
	string _thickEndStr;
	string _blockCountStr;
private:

	static inline void print_record(const Bed12Interval& what, CHRPOS start, CHRPOS end, string& outBuf);
};



#endif /* BED12INTERVAL_H_ */
