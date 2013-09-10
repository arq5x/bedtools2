#include "Bed6Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

Bed6Interval::Bed6Interval()
{

}

Bed6Interval::~Bed6Interval()
{

}

bool Bed6Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed3Interval::initFromFile(fileReader);

	fileReader->getField(3, _name);
	fileReader->getField(4, _score);
	char strandChar = 0;
	fileReader->getField(5, strandChar);
	setStrand(strandChar);

	return baseRetFlag;
}

void Bed6Interval::print(QuickString &outBuf) const
{
	Bed3Interval::print(outBuf);

	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(_score);
	outBuf.append('\t');
	outBuf.append(getStrandChar());
}

void Bed6Interval::print(QuickString &outBuf, int start, int end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(_score);
	outBuf.append('\t');
	outBuf.append(getStrandChar());
}

void Bed6Interval::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(_score);
	outBuf.append('\t');
	outBuf.append(getStrandChar());
}


void Bed6Interval::printNull(QuickString &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append("\t.\t-1\t.", 7);
}

