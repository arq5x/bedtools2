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

	fileReader->getField(3, nameRef());
	return baseRetFlag;
}

inline void print_record(const string& chrom_name, CHRPOS start, CHRPOS end, const string& name, string& buf) {
	append_cstr(buf, chrom_name);
	buf.push_back('\t');
	append_chrpos(buf, start);
	buf.push_back('\t');
	append_chrpos(buf, end);
	buf.push_back('\t');
	append_cstr(buf, name);
}

void BedGraphInterval::print(string &outBuf) const
{
	print_record(_chrName, _startPos, _endPos, nameRef(), outBuf);
}

void BedGraphInterval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(_chrName, start, end, nameRef(), outBuf);
}

void BedGraphInterval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	outBuf.append(nameRef());
}


void BedGraphInterval::printNull(string &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append("\t.", 2);
}

const string &BedGraphInterval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case 4:
		return nameRef();
		break;
	default:
		return Bed3Interval::getField(fieldNum);
		break;
	}
}

bool BedGraphInterval::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case 4:
		return true;
		break;
	default:
		return Bed3Interval::isNumericField(fieldNum);
		break;
	}
}

