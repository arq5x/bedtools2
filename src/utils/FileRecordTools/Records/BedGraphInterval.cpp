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

	fileReader->getField(3, _name);
	return baseRetFlag;
}

inline void print_record(const string& chrom_name, CHRPOS start, CHRPOS end, const string& name, string& buf) {
	const char* buffer = buffer_printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS"\t%s", chrom_name.c_str(), start, end, name.c_str());
	buf.append(buffer);
}

void BedGraphInterval::print(string &outBuf) const
{
	print_record(_chrName, _startPos, _endPos, _name, outBuf);
}

void BedGraphInterval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(_chrName, start, end, _name, outBuf);
}

void BedGraphInterval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	outBuf.append(_name);
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
		return _name;
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

