#include "BedPlusInterval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

BedPlusInterval::BedPlusInterval()
:  _numFixedFields(defaultNumFixedFields),
  _numPrintFields(0)
{
	_plusFields.setNumOffsetFields(defaultNumFixedFields);
}

void BedPlusInterval::setNumFixedFields(int numFields) {
	_numFixedFields = numFields;
	_plusFields.setNumOffsetFields(numFields);
}


bool BedPlusInterval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed3Interval::initFromFile(fileReader);
	if (fileReader->getCoordsOnly()) return baseRetFlag;

	if (_numFixedFields != defaultNumFixedFields) {
		if (!initNameScoreStrand(fileReader)) {
			fileReader->getField(3, nameRef());
			fileReader->getField(4, scoreRef());
			fileReader->getField(5, strandRef());
			adjustStrandVal();
		}
	}
	_plusFields.initFromFile(fileReader);
	return baseRetFlag;

}


bool BedPlusInterval::initFastLine(const ChrRef &chr, CHRPOS s, CHRPOS e, int fileIdx, const char *line, size_t lineLen, size_t tail, int numFields)
{
	if (numFields <= _numFixedFields) return false;
	if (_numFixedFields == defaultNumFixedFields) {
		setFastCoords(chr, s, e, fileIdx);
		_plusFields.setFast(line + tail + 1, lineLen - tail - 1, (size_t)(numFields - _numFixedFields));
		return true;
	}
	if (_numFixedFields != 6) return false;
	const char *p = line + tail + 1;
	const char *end = line + lineLen;
	const char *t = (const char *)memchr(p, '\t', end - p);
	if (t == NULL) return false;
	t = (const char *)memchr(t + 1, '\t', end - t - 1);
	if (t == NULL) return false;
	t = (const char *)memchr(t + 1, '\t', end - t - 1);
	setFastCoords(chr, s, e, fileIdx);
	if (t == NULL) {
		setFastNameScoreStrand(p, end - p);
		_plusFields.setFast(end, 0, 0);
	} else {
		setFastNameScoreStrand(p, t - p);
		_plusFields.setFast(t + 1, end - t - 1, (size_t)(numFields - 6));
	}
	return true;
}

void BedPlusInterval::clear() {
	Bed3Interval::clear();
	_plusFields.clear();
}

void BedPlusInterval::print(string &outBuf) const
{
	Bed3Interval::print(outBuf);
	outBuf.append("\t");
	printBed6PlusFields(outBuf);
	_plusFields.printFields(outBuf);
}

void BedPlusInterval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	printBed6PlusFields(outBuf);
	_plusFields.printFields(outBuf);
}

void BedPlusInterval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	printBed6PlusFields(outBuf);
	_plusFields.printFields(outBuf);
}


void BedPlusInterval::printNull(string &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	printBed6PlusNullFields(outBuf);
	for (int i=_numFixedFields; i < _numPrintFields; i++) {
		outBuf.append("\t.");
	}
}

const string &BedPlusInterval::getField(int fieldNum) const
{
	if (fieldNum > _numFixedFields) {
		return _plusFields.getField(fieldNum);
	} else if (fieldNum == 4 && _numFixedFields >=4) {
		return nameRef();
	} else if (fieldNum == 5 && _numFixedFields >=5) {
		return scoreRef();
	}
	else if (fieldNum == 6 && _numFixedFields >=6) {
		return strandRef();
	}
	return Bed3Interval::getField(fieldNum);
}

bool BedPlusInterval::isNumericField(int fieldNum) {

	//
	// TBD: There is no currently no good way to guarantee / enforce whether
	// fields after the 3rd are numeric, so for now we'll give the user the
	// benefit of the doubt on those.
	//
	if (fieldNum > defaultNumFixedFields) {
		return true;
	}
	return Bed3Interval::isNumericField(fieldNum);
}

void BedPlusInterval::printBed6PlusFields(string &outBuf) const {
	if (_numFixedFields != defaultNumFixedFields) {
		appendNameScoreStrand(outBuf);
		outBuf.append("\t");
	}
}

void BedPlusInterval::printBed6PlusNullFields(string &outBuf) const {
	if (_numFixedFields != defaultNumFixedFields) {
		outBuf.append("\t.\t.\t.");
	}

}
