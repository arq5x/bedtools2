#include "GffRecord.h"
#include "SingleLineDelimTextFileReader.h"
#include <cstring>

GffRecord::GffRecord() {

}

GffRecord::~GffRecord() {

}

void GffRecord::clear()
{
	Bed6Interval::clear();
	_source.clear();
	_frame.clear();
	_group.clear();
	_numFields = 8;
}

bool GffRecord::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	setFileIdx(fileReader->getFileIdx());
	fileReader->getField(0, _chrName);
	fileReader->getField(3, _startPosStr);
	_startPos = str2chrPos(_startPosStr);
	_startPos--; // VCF is one-based. Here we intentionally don't decrement the string version,
	//because we'll still want to output the one-based number in the print methods, even though
	//internally we decrement the integer to comply with the 0-based format common to other records.
	fileReader->getField(4, _endPosStr);
	//endPos is just the startPos plus the length of the variant
	_endPos = str2chrPos(_endPosStr);

	fileReader->getField(2, _name);
	fileReader->getField(1, _source);
	fileReader->getField(5, _score);

	//GFF allows a '.' for the strandChar, signifying it is not known.
	fileReader->getField(6, _strand);
	adjustStrandVal();

	fileReader->getField(7, _frame);
	_numFields = min(9, fileReader->getNumFields());
	if (_numFields == 9) {
		fileReader->getField(8, _group);
	}


	return true;
}

void GffRecord::print(string &outBuf) const
{
	ostringstream s;	
	s << _chrName;
	s << "\t";
	s << _source;
	s << "\t";
	s << _name;
	s << "\t";
	s << _startPosStr;
	s << "\t";
	s << _endPosStr;
	s << "\t";
	outBuf.append(s.str());
	printRemainingFields(outBuf);
}

void GffRecord::print(string &outBuf, CHRPOS start, CHRPOS end) const
{
	ostringstream s;	
	s << _chrName;
	s << "\t";
	s << _source;
	s << "\t";
	s << _name;
	s << "\t";
	s << start;
	s << "\t";
	s << end;
	s << "\t";
	outBuf.append(s.str());
	printRemainingFields(outBuf);
}

void GffRecord::print(string &outBuf, const string & start, const string & end) const
{
	outBuf.append(_chrName);
	outBuf.append("\t");
	outBuf.append(_source);
	outBuf.append("\t");
	outBuf.append(_name);
	outBuf.append("\t");
	outBuf.append(start);
	outBuf.append("\t");
	outBuf.append(end);
	outBuf.append("\t");

	printRemainingFields(outBuf);
}

void GffRecord::printRemainingFields(string &outBuf) const
{
	outBuf.append(_score);
	outBuf.append("\t");
	outBuf.append(_strand);
	outBuf.append("\t");
	outBuf.append(_frame);
	if (_numFields == 9) {
		outBuf.append("\t");
		outBuf.append(_group);
	}
}

void GffRecord::printNull(string &outBuf) const
{
	outBuf.append(".\t.\t.\t-1\t-1\t.\t.\t.", 17);
	if (_numFields > 8) {
		outBuf.append("\t.", 2);
	}
}

const string &GffRecord::getField(int fieldNum) const
{
	if (fieldNum == 9 && _numFields == 9) {
		return _group;
	}
	switch (fieldNum) {
	case 1:
		return _chrName;
		break;
	case 2:
		return _source;
		break;
	case 3:
		return _name;
		break;
	case 4:
		return _startPosStr;
		break;
	case 5:
		return _endPosStr;
		break;
	case 6:
		return _score;
		break;
	case 7:
		return _strand;
		break;
	case 8:
		return _frame;
		break;
	default:
		return Bed6Interval::getField(fieldNum);
		break;
	}
}

bool GffRecord::isNumericField(int fieldNum) {
	switch (fieldNum) {
	case 1:
		return false;
		break;
	case 2:
		return false;
		break;
	case 3:
		return false;
		break;
	case 4:
		return true;
		break;
	case 5:
		return true;
		break;
	case 6:
		return true;
		break;
	case 7:
		return false;
		break;
	case 8:
		return false;
		break;
	case 9:
		return false;
		break;
	default:
		return Bed6Interval::isNumericField(fieldNum);
		break;
	}

}


