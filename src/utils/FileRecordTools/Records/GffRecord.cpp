#include "GffRecord.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

GffRecord::GffRecord()
{

}

GffRecord::~GffRecord()
{

}

void GffRecord::clear()
{
	Bed6Interval::clear();
	_source.clear();
	_frame.clear();
	_group.clear();
	_numFields = 8;
}

bool GffRecord::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	fileReader->getField(0, _chrName);
	fileReader->getField(3, _startPosStr);
	_startPos = str2chrPos(_startPosStr);
	_startPos--; // VCF is one-based. Here we intentionally don't decrement the string version,
	//because we'll still want to output the one-based number in the print methods, even though
	//internally we decrement the integer to comply with the 0-based format common to other records.
	fileReader->getField(4, _endPosStr);
	//endPos is just the startPos plus the length of the variant
	_endPos = str2chrPos(_endPosStr);

	fileReader->getField(2, _name);
	fileReader->getField(1, _source);
	fileReader->getField(5, _score);

	//GFF allows a '.' for the strandChar, signifying it is not known.
	char strandChar = 0;
	fileReader->getField(6, strandChar);
	setStrand(strandChar);

	fileReader->getField(7, _frame);
	_numFields = fileReader->getNumFields();
	if (_numFields == 9) {
		fileReader->getField(8, _group);
	}


	return true;
}

void GffRecord::print(QuickString &outBuf) const
{
	outBuf.append(_chrName);
	outBuf.append('\t');
	outBuf.append(_source);
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(_startPosStr);
	outBuf.append('\t');
	outBuf.append(_endPosStr);
	outBuf.append('\t');

	printRemainingFields(outBuf);
}

void GffRecord::print(QuickString &outBuf, int start, int end) const
{
	outBuf.append(_chrName);
	outBuf.append('\t');
	outBuf.append(_source);
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(start +1);
	outBuf.append('\t');
	outBuf.append(end);
	outBuf.append('\t');

	printRemainingFields(outBuf);
}

void GffRecord::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	outBuf.append(_chrName);
	outBuf.append('\t');
	outBuf.append(_source);
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(start);
	outBuf.append('\t');
	outBuf.append(end);
	outBuf.append('\t');

	printRemainingFields(outBuf);
}

void GffRecord::printRemainingFields(QuickString &outBuf) const
{
	outBuf.append(_score);
	outBuf.append('\t');
	outBuf.append(getStrandChar());
	outBuf.append('\t');
	outBuf.append(_frame);
	if (_numFields == 9) {
		outBuf.append('\t');
		outBuf.append(_group);
	}
}

void GffRecord::printNull(QuickString &outBuf) const
{
	outBuf.append(".\t.\t.\t-1\t-1\t.\t.\t.", 17);
	if (_numFields > 8) {
		outBuf.append("\t.", 2);
	}
}

