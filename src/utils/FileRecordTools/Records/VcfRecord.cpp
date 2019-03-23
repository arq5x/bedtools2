/*
 * VcfRecord.cpp
 *
 *  Created on: May 1, 2013
 *      Author: nek3d
 */

#include "VcfRecord.h"
#include "SingleLineDelimTextFileReader.h"

bool VcfRecord::initFromFile(SingleLineDelimTextFileReader *fileReader)
{
	setFileIdx(fileReader->getFileIdx());
	fileReader->getField(0, _chrName);
	_chrId = fileReader->getCurrChromdId();
	fileReader->getField(1, _startPosStr);
	_startPos = str2chrPos(_startPosStr);
	_startPos--; // VCF is one-based. Here we intentionally don't decrement the string version,
	//because we'll still want to output the one-based number in the print methods, even though
	//internally we decrement the integer to comply with the 0-based format common to other records.
	fileReader->getField(4, _varAlt);
	fileReader->getField(3, _varRef);
	if (_varAlt[0] == '<') {
		// this is a structural variant. Need to parse the tags to find the endpos,
		// UNLESS it's an insertion.
		if (!(_varAlt[1] == 'I' && _varAlt[2] == 'N' && _varAlt[3] == 'S')) {
			_endPos = _startPos + fileReader->getVcfSVlen();
		} else {
			//for insertions, treat as zero-length records
			_endPos = _startPos;
		}
	} else {
		//endPos is just the startPos plus the length of the variant
		_endPos = _startPos + _varRef.size();
	}
	int2str(_endPos, _endPosStr);
	fileReader->getField(2, _name);
	fileReader->getField(5, _score);

	_plusFields.setNumOffsetFields(numFixedFields);
	return _plusFields.initFromFile(fileReader);
}

void VcfRecord::clear()
{
	BedPlusInterval::clear();
	_varRef.clear();
	_varAlt.clear();
}

void VcfRecord::print(string &outBuf) const {
	outBuf.append(_chrName);
	outBuf.append("\t");
	outBuf.append(_startPosStr);
	printOtherFields(outBuf);
}

void VcfRecord::print(string &outBuf, CHRPOS start, CHRPOS end) const {
	outBuf.append(_chrName);
	outBuf.append("\t");
	outBuf.append(_startPosStr);
	printOtherFields(outBuf);
}

void VcfRecord::print(string &outBuf, const string & start, const string & end) const {
	outBuf.append(_chrName);
	outBuf.append("\t");
	outBuf.append(_startPosStr);
	printOtherFields(outBuf);

}

void VcfRecord::printNull(string &outBuf) const {
	outBuf.append(".\t-1\t.");
	for (int i=3; i < _numPrintFields; i++) {
		outBuf.append("\t.");
	}
}

void VcfRecord::printOtherFields(string &outBuf) const {
	outBuf.append("\t");
	outBuf.append(_name);
	outBuf.append("\t");
	outBuf.append(_varRef);
	outBuf.append("\t");
	outBuf.append(_varAlt);
	outBuf.append("\t");
	outBuf.append(_score);
	outBuf.append("\t");
	_plusFields.printFields(outBuf);
}

const string &VcfRecord::getField(int fieldNum) const
{
	//a request for any of the first six fields will retrieve
	//chrom, start, name, varRef, varAlt, score,  in that order.

	switch (fieldNum) {
	case 1:
		return _chrName;
		break;
	case 2:
		return _startPosStr;
		break;
	case 3:
		return _name;
	case 4:
		return _varRef;
		break;
	case 5:
		return _varAlt;
		break;
	case 6:
		return _score;
	default:
		return BedPlusInterval::getField(fieldNum);
		break;
	}
}

bool VcfRecord::isNumericField(int fieldNum) {

	switch (fieldNum) {
	case 1:
		return false; //_chrName;
		break;
	case 2:
		return true; //_startPosStr;
		break;
	case 3:
		return false; //_name;
	case 4:
		return false; //_varRef;
		break;
	case 5:
		return false; //_varAlt;
		break;
	case 6:
		return true; //_score;
	default:
		return BedPlusInterval::isNumericField(fieldNum);
		break;
	}
}
