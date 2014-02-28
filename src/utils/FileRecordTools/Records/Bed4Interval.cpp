#include "Bed4Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

Bed4Interval::Bed4Interval()
{

}

Bed4Interval::~Bed4Interval()
{

}

bool Bed4Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed3Interval::initFromFile(fileReader);

	fileReader->getField(3, _name);
	return baseRetFlag;
}

void Bed4Interval::print(QuickString &outBuf) const
{
	Bed3Interval::print(outBuf);
	outBuf.append('\t');
	outBuf.append(_name);
}

void Bed4Interval::print(QuickString &outBuf, int start, int end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_name);
}

void Bed4Interval::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_name);
}

void Bed4Interval::printNull(QuickString &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append('\t');
	outBuf.append("\t.", 2);
}

const QuickString &Bed4Interval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case 4:
		return _name;
		break;
	default:
		return Bed3Interval::getField(fieldNum);
		break;
	}
}

bool Bed4Interval::isNumericField(int fieldNum) {
	return (fieldNum == 4 ? false : Bed3Interval::isNumericField(fieldNum));
}


