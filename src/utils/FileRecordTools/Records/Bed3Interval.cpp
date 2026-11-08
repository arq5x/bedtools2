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
	{
		const char *line = fileReader->getLineData();
		const int *dp = fileReader->getDelimPositions();
		_chrName.assign(line + dp[0] + 1, dp[1] - dp[0] - 1);
	}
	string s1, s2;
	bool canonStart = fileReader->getPosField(1, _startPos, s1);
	bool canonEnd = fileReader->getPosField(2, _endPos, s2);
	if (canonStart && canonEnd) {
		_posMode = POS_PLAIN;
		_origStart = _startPos;
		_origEnd = _endPos;
	} else {
		// not the plain decimal form (leading zeros, exponent, ...): keep the original text
		if (canonStart) fileReader->getField(1, s1);
		if (canonEnd) fileReader->getField(2, s2);
		PosStrBox *b = _posBox.make();
		b->start.swap(s1);
		b->end.swap(s2);
		_posMode = POS_STRINGS;
	}
	_zeroLength = false;
	return true;
}

bool Bed3Interval::initFastLine(const ChrRef &chr, CHRPOS s, CHRPOS e, int fileIdx, const char *, size_t, size_t, int)
{
	setFastCoords(chr, s, e, fileIdx);
	return true;
}

inline void print_record(const string& name, CHRPOS start, CHRPOS end, string& buf) {
	append_cstr(buf, name);
	buf.push_back('\t');
	append_chrpos(buf, start);
	buf.push_back('\t');
	append_chrpos(buf, end);
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
		return _chrName.str();
		break;
	case 2:
		return startStr();
		break;
	case 3:
		return endStr();
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
