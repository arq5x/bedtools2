#include "BedPlusInterval.h"
#include "SingleLineDelimTextFileReader.h"

BedPlusInterval::BedPlusInterval()
:  _numPrintFields(0)
{
	_plusFields.setNumOffsetFields(numFixedFields);
}


bool BedPlusInterval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	return (Bed3Interval::initFromFile(fileReader) && _plusFields.initFromFile(fileReader));
}


void BedPlusInterval::clear() {
	Bed3Interval::clear();
	_plusFields.clear();
}

void BedPlusInterval::print(QuickString &outBuf) const
{
	Bed3Interval::print(outBuf);
	_plusFields.printFields(outBuf);
}

void BedPlusInterval::print(QuickString &outBuf, int start, int end) const
{
	Bed3Interval::print(outBuf, start, end);
	_plusFields.printFields(outBuf);
}

void BedPlusInterval::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	Bed3Interval::print(outBuf, start, end);
	_plusFields.printFields(outBuf);
}


void BedPlusInterval::printNull(QuickString &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	for (int i=numFixedFields; i < _numPrintFields; i++) {
		outBuf.append("\t.");
	}
}

const QuickString &BedPlusInterval::getField(int fieldNum) const
{
	if (fieldNum > numFixedFields) {
		return _plusFields.getField(fieldNum);
	}
	return Bed3Interval::getField(fieldNum);
}

bool BedPlusInterval::isNumericField(int fieldNum) {

	//
	// TBD: There is no currently no good way to guarantee / enforce whether
	// fields after the 3rd are numeric, so for now we'll give the user the
	// benefit of the doubt on those.
	//
	if (fieldNum > numFixedFields) {
		return true;
	}
	return Bed3Interval::isNumericField(fieldNum);
}

