#include "Bed3Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstdio>
#include <cstdlib>

Bed3Interval::Bed3Interval()
{

}

Bed3Interval::~Bed3Interval()
{

}

bool Bed3Interval::initFromFile(FileReader *fileReader)
{
	SingleLineDelimTextFileReader *sldFileReader = static_cast<SingleLineDelimTextFileReader*>(fileReader);
	return initFromFile(sldFileReader);
}


bool Bed3Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	setFileIdx(fileReader->getFileIdx());
	fileReader->getField(0, _chrName);
	fileReader->getField(1, _startPosStr);
	fileReader->getField(2, _endPosStr);
	_startPos = str2chrPos(_startPosStr);
	_endPos = str2chrPos(_endPosStr);
	_zeroLength = false;
	return true;
}

inline void print_record(const string& name, CHRPOS start, CHRPOS end, string& buf) {
	const char* buffer = buffer_printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS, name.c_str(), start, end);
	buf.append(buffer);
}

void Bed3Interval::print(string &outBuf) const{
	print_record(_chrName, _startPos, _endPos, outBuf);
}

void Bed3Interval::print(string &outBuf, CHRPOS start, CHRPOS end) const{
	print_record(_chrName, start, end, outBuf);
}

void Bed3Interval::print(string &outBuf, const string & start, const string & end) const
{
	outBuf.append(_chrName);
	outBuf.append("\t");
	outBuf.append(start);
	outBuf.append("\t");
	outBuf.append(end);
}

void Bed3Interval::printNull(string &outBuf) const {
	outBuf.append(".\t-1\t-1", 7);
}

const string &Bed3Interval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case 1:
		return _chrName;
		break;
	case 2:
		return _startPosStr;
		break;
	case 3:
		return _endPosStr;
		break;
	default:
		return Record::getField(fieldNum);
		break;
	}
}

bool Bed3Interval::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case 1:
		return false; //chrom
		break;
	case 2:
		return true; //startPos
		break;
	case 3:
		return true; //endPos
		break;
	default:
	    cerr << endl << "*****" << endl
	         << "*****ERROR: requested invalid column " << fieldNum << ". Exiting." << endl
	          << endl << "*****" << endl;
	    exit(1);
	    break;
	}
}
