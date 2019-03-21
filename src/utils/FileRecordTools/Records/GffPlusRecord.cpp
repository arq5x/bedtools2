#include "GffPlusRecord.h"
#include "SingleLineDelimTextFileReader.h"

GffPlusRecord::GffPlusRecord()
: _numPrintFields(0)
{
}

GffPlusRecord::~GffPlusRecord() {

}

bool GffPlusRecord::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	if (!GffRecord::initFromFile(fileReader)) {
		return false;
	}
	_plusFields.setNumOffsetFields(GffRecord::getNumFields());
	return _plusFields.initFromFile(fileReader);
}


void GffPlusRecord::clear() {
	GffRecord::clear();
	_plusFields.clear();
}

void GffPlusRecord::print(string &outBuf) const
{
	GffRecord::print(outBuf);
	_plusFields.printFields(outBuf);
}

void GffPlusRecord::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	GffRecord::print(outBuf, start, end);
	_plusFields.printFields(outBuf);
}

void GffPlusRecord::print(string &outBuf, const string & start, const string & end) const
{
	GffRecord::print(outBuf, start, end);
	_plusFields.printFields(outBuf);
}


void GffPlusRecord::printNull(string &outBuf) const
{
	GffRecord::printNull(outBuf);
	for (int i=_numFields; i < _numPrintFields; i++) {
		outBuf.append("\t.");
	}
}

const string &GffPlusRecord::getField(int fieldNum) const
{
	if (fieldNum > _numFields) {
		return _plusFields.getField(fieldNum);
	}
	return GffRecord::getField(fieldNum);
}

bool GffPlusRecord::isNumericField(int fieldNum) {

	if (fieldNum < 9) {
		return GffRecord::isNumericField(fieldNum);
	}
	return true;
}

