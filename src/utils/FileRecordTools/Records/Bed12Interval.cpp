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
	 _blockCount = (int)str2chrPos(_blockCountStr);
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

 inline void Bed12Interval::print_record(const Bed12Interval& what, CHRPOS start, CHRPOS end, string& outBuf) {
	const char* buffer = buffer_printf(
			"%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS"\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
			what._chrName.c_str(), 
			start, 
			end, 
			what._name.c_str(), 
			what._score.c_str(), 
			what._strand.c_str(),
			what._thickStartStr.c_str(),
			what._thickEndStr.c_str(),
			what._itemRGB.c_str(),
			what._blockCountStr.c_str(),
			what._blockSizes.c_str(),
			what._blockStarts.c_str()
	);

	outBuf.append(buffer);
}

void Bed12Interval::print(string &outBuf) const
{
	print_record(*this, _startPos, _endPos, outBuf);
}

void Bed12Interval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(*this, start, end,outBuf);
}

void Bed12Interval::print(string &outBuf, const string & start, const string & end) const
{
	Bed6Interval::print(outBuf, start, end);
	outBuf.append("\t");
	outBuf.append(_thickStartStr);
	outBuf.append("\t");
	outBuf.append(_thickEndStr);
	outBuf.append("\t");
	outBuf.append(_itemRGB);
	outBuf.append("\t");
	outBuf.append(_blockCountStr);
	outBuf.append("\t");
	outBuf.append(_blockSizes);
	outBuf.append("\t");
	outBuf.append(_blockStarts);
}


void Bed12Interval::printNull(string &outBuf) const
{
	Bed6Interval::printNull(outBuf);

	outBuf.append("\t.\t.\t.\t.\t.\t.", 12);
}

const string &Bed12Interval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case 7:
		return _thickStartStr;
		break;
	case 8:
		return _thickEndStr;
		break;
	case 9:
		return _itemRGB;
		break;
	case 10:
		return _blockCountStr;
		break;
	case 11:
		return _blockSizes;
		break;
	case 12:
		return _blockStarts;
		break;
	default:
		return Bed6Interval::getField(fieldNum);
		break;
	}
}

bool Bed12Interval::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case 7:
		return true;
		break;
	case 8:
		return true;
		break;
	case 9:
		return false;
		break;
	case 10:
		return true;
		break;
	case 11:
		return false;
		break;
	case 12:
		return false;
		break;
	default:
		return Bed6Interval::isNumericField(fieldNum);
		break;
	}
}

CHRPOS Bed12Interval::getLength(bool obeySplits) const {
	//only bed12 and BAM need to check splits
	if (!obeySplits || _blockCount <=0) {
		return _endPos - _startPos;
	} else {
		vector<CHRPOS> vBlockSizes;
		Tokenize(_blockSizes, vBlockSizes, ',');
	    return accumulate(vBlockSizes.begin(), vBlockSizes.end(), 0);
	}
}
