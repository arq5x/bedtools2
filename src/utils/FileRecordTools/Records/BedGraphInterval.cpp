#include "BedGraphInterval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

BedGraphInterval::BedGraphInterval()
{

}

BedGraphInterval::~BedGraphInterval()
{

}

bool BedGraphInterval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed3Interval::initFromFile(fileReader);

	fileReader->getField(3, _name);
	return baseRetFlag;
}

void BedGraphInterval::print(QuickString &outBuf) const
{
	Bed3Interval::print(outBuf);
	outBuf.append('\t');
	outBuf.append(_name);
}

void BedGraphInterval::print(QuickString &outBuf, int start, int end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_name);
}

void BedGraphInterval::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_name);
}


void BedGraphInterval::printNull(QuickString &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append("\t.", 2);
}

