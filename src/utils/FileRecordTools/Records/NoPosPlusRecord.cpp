#include "NoPosPlusRecord.h"
#include "SingleLineDelimTextFileReader.h"

NoPosPlusRecord::NoPosPlusRecord()
{
	_plusFields.setNumOffsetFields(defaultNumFixedFields);
}

NoPosPlusRecord::~NoPosPlusRecord() {

}

bool NoPosPlusRecord::initFromFile(FileReader *fileReader)
{
	SingleLineDelimTextFileReader *sldFileReader = static_cast<SingleLineDelimTextFileReader*>(fileReader);
	return initFromFile(sldFileReader);
}

bool NoPosPlusRecord::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	return _plusFields.initFromFile(fileReader);
}


void NoPosPlusRecord::clear() {
	_plusFields.clear();
}

void NoPosPlusRecord::print(string &outBuf) const
{
	_plusFields.printFields(outBuf);
}

const string &NoPosPlusRecord::getField(int fieldNum) const
{
	return _plusFields.getField(fieldNum);
}

bool NoPosPlusRecord::isNumericField(int fieldNum) {
	return true;
}

