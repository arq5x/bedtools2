#include "BedPlusInterval.h"
#include "SingleLineDelimTextFileReader.h"

BedPlusInterval::BedPlusInterval()
: _numPrintFields(0)
{
}

BedPlusInterval::~BedPlusInterval()
{
	for (int i=0; i < (int)_otherIdxs.size(); i++) {
		delete _otherIdxs[i];
	}
}

const BedPlusInterval &BedPlusInterval::operator=(const BedPlusInterval &other) {
	Bed6Interval::operator=(other);

	int otherSize = other._otherIdxs.size();
	int mySize = _otherIdxs.size();

	_numPrintFields = other._numPrintFields;
	int numMatchingFields = min(mySize, otherSize);
	for (int i=0; i < numMatchingFields; i++) {
		(*(_otherIdxs[i])) = (*(other._otherIdxs[i]));
	}
	if (mySize < otherSize) {
		for (int i = mySize; i < otherSize; i++) {
			QuickString *pqs = new QuickString(*(other._otherIdxs[i]));
			_otherIdxs.push_back(pqs);
		}
	} else if (mySize > otherSize) {
		for (int i= otherSize; i < mySize; i++) {
			delete _otherIdxs[i];
		}
		_otherIdxs.resize(otherSize);
	}
	return *this;
}

bool BedPlusInterval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	return (Bed6Interval::initFromFile(fileReader) && initOtherFieldsFromFile(fileReader));
}

bool BedPlusInterval::initOtherFieldsFromFile(SingleLineDelimTextFileReader *fileReader)
{
	int numFields = fileReader->getNumFields() - startOtherIdx;

	if ((int)_otherIdxs.size() != numFields) {
		if ((int)_otherIdxs.size() > 0) {
			return false; //file had a number of fields not matching what was expected.
		}
		for (int i=0; i < numFields; i++) {
			_otherIdxs.push_back(new QuickString());
		}
	}

	for (int i=0; i < numFields; i++) {
		fileReader->getField(i + startOtherIdx, (*(_otherIdxs[i])));
	}
	return true;
}

void BedPlusInterval::clear() {
	Bed6Interval::clear();
	_numPrintFields = 0;
	for (int i=0; i < (int)_otherIdxs.size(); i++) {
		_otherIdxs[i]->clear();
	}
}

void BedPlusInterval::print(QuickString &outBuf) const
{
	Bed6Interval::print(outBuf);

	for (int i=0; i < (int)_otherIdxs.size(); i++) {
		outBuf.append('\t');
		outBuf.append(*(_otherIdxs[i]));
	}
}

void BedPlusInterval::print(QuickString &outBuf, int start, int end) const
{
	Bed6Interval::print(outBuf, start, end);
	for (int i=0; i < (int)_otherIdxs.size(); i++) {
		outBuf.append('\t');
		outBuf.append(*(_otherIdxs[i]));
	}
}

void BedPlusInterval::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const
{
	Bed6Interval::print(outBuf, start, end);
	for (int i=0; i < (int)_otherIdxs.size(); i++) {
		outBuf.append('\t');
		outBuf.append(*(_otherIdxs[i]));
	}
}


void BedPlusInterval::printNull(QuickString &outBuf) const
{
	Bed6Interval::printNull(outBuf);
	for (int i=startOtherIdx; i < _numPrintFields; i++) {
		outBuf.append("\t.");
	}
}

const QuickString &BedPlusInterval::getField(int fieldNum) const
{
	//a request for any of the first six fields will retrieve
	//chrom, start, end, name, score, and strand, in that order.
	//A request for field 6+ will go to the otherIdxs.
	if (fieldNum > startOtherIdx && fieldNum <= startOtherIdx + (int)_otherIdxs.size()) {
		return (*(_otherIdxs[fieldNum - startOtherIdx - 1]));
	}
	return Bed6Interval::getField(fieldNum);
}

bool BedPlusInterval::isNumericField(int fieldNum) {

	//
	// TBD: There is no currently no good way to guarantee / enforce whether
	// fields after the 6th are numeric, so for now we'll give the user the
	// benefit of the doubt on those.
	//
	if (fieldNum > startOtherIdx) {
		return true;
	} else {
		return Bed6Interval::isNumericField(fieldNum);
	}
}

