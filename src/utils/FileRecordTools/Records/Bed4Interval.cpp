#include "Bed4Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

Bed4Interval::Bed4Interval()
{

}

Bed4Interval::~Bed4Interval()
{

}

bool Bed4Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
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

void Bed4Interval::print(string &outBuf) const
{
	print_record(_chrName, _startPos, _endPos, nameRef(), outBuf);
}

void Bed4Interval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(_chrName, start, end, nameRef(), outBuf);
}

void Bed4Interval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	outBuf.append(nameRef());
}

void Bed4Interval::printNull(string &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append("\t");
	outBuf.append(".");
}

const string &Bed4Interval::getField(int fieldNum) const
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

bool Bed4Interval::isNumericField(int fieldNum) {
	return (fieldNum == 4 ? false : Bed3Interval::isNumericField(fieldNum));
}


