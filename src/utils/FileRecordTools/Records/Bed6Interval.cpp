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
	if (fileReader->getCoordsOnly()) return baseRetFlag;

	if (!initNameScoreStrand(fileReader)) {
		fileReader->getField(3, nameRef());
		fileReader->getField(4, scoreRef());
		fileReader->getField(5, strandRef());
		adjustStrandVal();
	}
	return baseRetFlag;
}

bool Bed6Interval::initFastLine(const ChrRef &chr, CHRPOS s, CHRPOS e, int fileIdx, const char *line, size_t lineLen, size_t tail, int numFields)
{
	if (numFields != 6 || tail >= lineLen) return false;
	setFastCoords(chr, s, e, fileIdx);
	setFastNameScoreStrand(line + tail + 1, lineLen - tail - 1);
	return true;
}

inline void print_record(const Record &rec, const string& chrom_name, CHRPOS start, CHRPOS end, string& buf) {
	append_cstr(buf, chrom_name);
	buf.push_back('\t');
	append_chrpos(buf, start);
	buf.push_back('\t');
	append_chrpos(buf, end);
	buf.push_back('\t');
	rec.appendNameScoreStrand(buf);
}

void Bed6Interval::print(string &outBuf) const
{
	print_record(*this, _chrName, _startPos, _endPos, outBuf);
}

void Bed6Interval::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	print_record(*this, _chrName, start, end, outBuf);
}

void Bed6Interval::print(string &outBuf, const string & start, const string & end) const
{
	Bed3Interval::print(outBuf, start, end);
	outBuf.push_back('\t');
	appendNameScoreStrand(outBuf);
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
		return nameRef();
		break;
	case 5:
		return scoreRef();
		break;
	case 6:
		return strandRef();
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
