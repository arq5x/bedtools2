#include "Bed5Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

Bed5Interval::Bed5Interval()
{

}

Bed5Interval::~Bed5Interval()
{

}

bool Bed5Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed3Interval::initFromFile(fileReader);

	fileReader->getField(3, _name);
	fileReader->getField(4, _score);
	return baseRetFlag;
}

void Bed5Interval::print(QuickString &outBuf) const
{
	Bed3Interval::print(outBuf);
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(_score);
}

void Bed5Interval::print(QuickString &outBuf, int start, int end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(_score);
}

void Bed5Interval::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(_score);
}


void Bed5Interval::printNull(QuickString &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append("\t.\t-1", 5);
}

const QuickString &Bed5Interval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case 4:
		return _name;
		break;
	case 5:
		return _score;
		break;

	default:
		return Bed3Interval::getField(fieldNum);
		break;
	}
}

bool Bed5Interval::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case 4:
		return false;
		break;
	case 5:
		return true;
		break;
	default:
		return Bed3Interval::isNumericField(fieldNum);
	}
}
