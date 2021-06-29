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

	fileReader->getField(3, _name);
	return baseRetFlag;
}

inline void print_record(const string& chrom_name, CHRPOS start, CHRPOS end, const string& name, string& buf) {
	const char* buffer = buffer_printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS"\t%s", chrom_name.c_str(), start, end, name.c_str());
	buf.append(buffer);
}

void Bed4Interval::print(string &outBuf) const
{
	print_record(_chrName, _startPos, _endPos, _name, outBuf);
}

void Bed4Interval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(_chrName, start, end, _name, outBuf);
}

void Bed4Interval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	outBuf.append(_name);
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
		return _name;
		break;
	default:
		return Bed3Interval::getField(fieldNum);
		break;
	}
}

bool Bed4Interval::isNumericField(int fieldNum) {
	return (fieldNum == 4 ? false : Bed3Interval::isNumericField(fieldNum));
}


