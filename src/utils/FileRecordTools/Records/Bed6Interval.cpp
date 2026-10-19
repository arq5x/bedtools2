#include "Bed6Interval.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

Bed6Interval::Bed6Interval()
{

}

Bed6Interval::~Bed6Interval()
{

}

bool Bed6Interval::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	bool baseRetFlag = Bed3Interval::initFromFile(fileReader);

	fileReader->getField(3, _name);
	fileReader->getField(4, _score);
	fileReader->getField(5, _strand);
	adjustStrandVal();
	return baseRetFlag;
}

inline void print_record(const string& chrom_name, CHRPOS start, CHRPOS end, const string& name, const string& score, const string& strand, string& buf) {
	append_cstr(buf, chrom_name);
	append_coords(buf, start, end);
	buf.push_back('\t');
	append_cstr(buf, name);
	buf.push_back('\t');
	append_cstr(buf, score);
	buf.push_back('\t');
	append_cstr(buf, strand);
}

void Bed6Interval::print(string &outBuf) const
{
	print_record(_chrName, _startPos, _endPos, _name, _score, _strand, outBuf);
}

void Bed6Interval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(_chrName, start, end, _name, _score, _strand, outBuf);
}

void Bed6Interval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.append("\t");
	outBuf.append(_name);
	outBuf.append("\t");
	outBuf.append(_score);
	outBuf.append("\t");
	outBuf.append(_strand);
}


void Bed6Interval::printNull(string &outBuf) const
{
	Bed3Interval::printNull(outBuf);
	outBuf.append("\t.\t-1\t.", 7);
}

const string &Bed6Interval::getField(int fieldNum) const
{
	switch (fieldNum) {
	case 4:
		return _name;
		break;
	case 5:
		return _score;
		break;
	case 6:
		return _strand;
		break;
	default:
		return Bed3Interval::getField(fieldNum);
		break;
	}
}

bool Bed6Interval::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case 4:
		return false;
		break;
	case 5:
		return true;
		break;
	case 6:
		return false;
		break;
	default:
		return Bed3Interval::isNumericField(fieldNum);
		break;
	}
}
