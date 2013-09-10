#include "Bed12Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstdlib>

Bed12Interval::Bed12Interval()
:_thickStart(-1),
 _thickEnd(-1),
 _blockCount(-1)

{
}

Bed12Interval::~Bed12Interval()
{
}
const Bed12Interval &Bed12Interval::operator=(const Bed12Interval &other) {
	Bed6Interval::operator=(other);

	_thickStart = other._thickStart;
	_thickEnd = other._thickEnd;
	_itemRGB = other._itemRGB;
	_blockCount = other._blockCount;
	_blockSizes = other._blockSizes;
	_blockStarts = other._blockStarts;

	return *this;
}

bool Bed12Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed6Interval::initFromFile(fileReader);

	 fileReader->getField(6, _thickStartStr);
	 fileReader->getField(7, _thickEndStr);
	 fileReader->getField(8, _itemRGB);
	 fileReader->getField(9, _blockCountStr);
	 fileReader->getField(10, _blockSizes);
	 fileReader->getField(11, _blockStarts);

	 _thickStart = str2chrPos(_thickStartStr);
	 _thickEnd = str2chrPos(_thickEndStr);
	 _blockCount = str2chrPos(_blockCountStr);
	return baseRetFlag;
}

void Bed12Interval::clear() {
	Bed6Interval::clear();

	 _thickStart = -1;
	 _thickEnd = -1;
	 _itemRGB.clear();
	 _blockCount = -1;
	 _blockSizes.clear();
	 _blockStarts.clear();
	 _thickStartStr.clear();
	 _thickEndStr.clear();
	 _blockCountStr.clear();

}

void Bed12Interval::print(QuickString &outBuf) const
{
	Bed6Interval::print(outBuf);

	outBuf.append('\t');
	outBuf.append(_thickStartStr);
	outBuf.append('\t');
	outBuf.append(_thickEndStr);
	outBuf.append('\t');
	outBuf.append(_itemRGB);
	outBuf.append('\t');
	outBuf.append(_blockCountStr);
	outBuf.append('\t');
	outBuf.append(_blockSizes);
	outBuf.append('\t');
	outBuf.append(_blockStarts);
}

void Bed12Interval::print(QuickString &outBuf, int start, int end) const
{
	Bed6Interval::print(outBuf, start, end);

	outBuf.append('\t');
	outBuf.append(_thickStartStr);
	outBuf.append('\t');
	outBuf.append(_thickEndStr);
	outBuf.append('\t');
	outBuf.append(_itemRGB);
	outBuf.append('\t');
	outBuf.append(_blockCountStr);
	outBuf.append('\t');
	outBuf.append(_blockSizes);
	outBuf.append('\t');
	outBuf.append(_blockStarts);
}

void Bed12Interval::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	Bed6Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_thickStartStr);
	outBuf.append('\t');
	outBuf.append(_thickEndStr);
	outBuf.append('\t');
	outBuf.append(_itemRGB);
	outBuf.append('\t');
	outBuf.append(_blockCountStr);
	outBuf.append('\t');
	outBuf.append(_blockSizes);
	outBuf.append('\t');
	outBuf.append(_blockStarts);
}


void Bed12Interval::printNull(QuickString &outBuf) const
{
	Bed6Interval::printNull(outBuf);

	outBuf.append("\t.\t.\t.\t.\t.\t.", 12);
}


