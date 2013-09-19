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

