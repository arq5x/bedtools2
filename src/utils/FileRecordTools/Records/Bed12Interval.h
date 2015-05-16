/*
 * Bed12Interval.h
 *
 *  Created on: Nov 28, 2012
 *      Author: nek3d
 */

#ifndef BED12INTERVAL_H_
#define BED12INTERVAL_H_

#include "Bed6Interval.h"

class SingleLineDelimTextFileReader;


class Bed12Interval : public Bed6Interval {
public:
	friend class FreeList<Bed12Interval>;

	Bed12Interval();
	virtual bool initFromFile(SingleLineDelimTextFileReader *);
	virtual void clear();
	virtual void print(QuickString &outBuf) const;
	virtual void print(QuickString &outBuf, int start, int end) const;
	virtual void print(QuickString &outBuf, const QuickString & start, const QuickString & end) const;
	virtual void printNull(QuickString &outBuf) const;
	virtual FileRecordTypeChecker::RECORD_TYPE getType() const { return FileRecordTypeChecker::BED12_RECORD_TYPE; }

	const Bed12Interval &operator=(const Bed12Interval &other);

	virtual int getThickStart() const { return _thickStart; }
	virtual void setThickStart(int thickStart)  { _thickStart = thickStart; }

	virtual int getThickEnd() const { return _thickEnd; }
	virtual void setThickEnd(int thickEnd) { _thickEnd = thickEnd; }

	virtual const QuickString & getItemRGB() const { return _itemRGB; }
	virtual void setItemRGB(const QuickString & rgb) { _itemRGB = rgb; }
	virtual void setItemRGB(const string & rgb) { _itemRGB = rgb; }
	virtual void setItemRGB(const char *rgb) { _itemRGB = rgb; }

	virtual int getBlockCount() const { return _blockCount; }
	virtual void setBlockCount(int blockCount) { _blockCount = blockCount; }

	virtual const QuickString & getBlockSizes() const { return _blockSizes; }
	virtual void setBlockSizes(const QuickString & blockSizes) { _blockSizes = blockSizes; }
	virtual void setBlockSizes(const string & blockSizes) { _blockSizes = blockSizes; }
	virtual void setBlockSizes(const char *blockSizes) { _blockSizes = blockSizes; }

	virtual const QuickString & getBlockStarts() const { return _blockStarts; }
	virtual void setBlockStarts(const QuickString & blockStarts) { _blockStarts = blockStarts; }
	virtual void setBlockStarts(const string & blockStarts) { _blockStarts = blockStarts; }
	virtual void setBlockStarts(const char *blockStarts) { _blockStarts = blockStarts; }

	virtual const QuickString &getField(int fieldNum) const;
	virtual int getNumFields() const  { return 12; }
	static bool isNumericField(int fieldNum);
	int getLength(bool obeySplits) const;


protected:
	virtual ~Bed12Interval();

	int _thickStart;
	int _thickEnd;
	QuickString _itemRGB;
	int _blockCount;
	QuickString _blockSizes;
	QuickString _blockStarts;

	//store int values as their originial strings for faster outputting.
	QuickString _thickStartStr;
	QuickString _thickEndStr;
	QuickString _blockCountStr;
};



#endif /* BED12INTERVAL_H_ */
