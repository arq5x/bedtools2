#include "Bed5Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

Bed5Interval::Bed5Interval()
{

}

Bed5Interval::~Bed5Interval()
{

}

bool Bed5Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed3Interval::initFromFile(fileReader);

	fileReader->getField(3, nameRef());
	fileReader->getField(4, scoreRef());
	return baseRetFlag;
}

inline void print_record(const string& chrom_name, CHRPOS start, CHRPOS end, const string& name, const string& score, string& buf) {
	append_cstr(buf, chrom_name);
	buf.push_back('\t');
	append_chrpos(buf, start);
	buf.push_back('\t');
	append_chrpos(buf, end);
	buf.push_back('\t');
	append_cstr(buf, name);
	buf.push_back('\t');
	append_cstr(buf, score);
}

void Bed5Interval::print(string &outBuf) const
{
	print_record(_chrName, _startPos, _endPos, nameRef(), scoreRef(), outBuf);
}

void Bed5Interval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(_chrName, start, end, nameRef(), scoreRef(), outBuf);
}

void Bed5Interval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	outBuf.append(nameRef());
	outBuf.append("\t");
	outBuf.append(scoreRef());
}


void Bed5Interval::printNull(string &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append("\t.\t-1", 5);
}

const string &Bed5Interval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case 4:
		return nameRef();
		break;
	case 5:
		return scoreRef();
		break;

	default:
		return Bed3Interval::getField(fieldNum);
		break;
	}
}

bool Bed5Interval::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case 4:
		return false;
		break;
	case 5:
		return true;
		break;
	default:
		return Bed3Interval::isNumericField(fieldNum);
	}
}
