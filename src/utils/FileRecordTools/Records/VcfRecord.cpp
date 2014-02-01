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
	fileReader->getField(0, _chrName);
	_chrId = fileReader->getCurrChromdId();
	fileReader->getField(1, _startPosStr);
	_startPos = str2chrPos(_startPosStr);
	_startPos--; // VCF is one-based. Here we intentionally don't decrement the string version,
	//because we'll still want to output the one-based number in the print methods, even though
	//internally we decrement the integer to comply with the 0-based format common to other records.
	fileReader->getField(3, _varAlt);
	//endPos is just the startPos plus the length of the variant
	_endPos = _startPos + _varAlt.size();
	int2str(_endPos, _endPosStr);

	fileReader->getField(2, _name);
	fileReader->getField(4, _varRef);
	fileReader->getField(5, _score);

	return initOtherFieldsFromFile(fileReader);
}

void VcfRecord::clear()
{
	BedPlusInterval::clear();
	_varAlt.clear();
	_varRef.clear();
}

void VcfRecord::print(QuickString &outBuf) const {
	outBuf.append(_chrName);
	outBuf.append('\t');
	outBuf.append(_startPosStr);
	printOtherFields(outBuf);
}

void VcfRecord::print(QuickString &outBuf, int start, int end) const {
	outBuf.append(_chrName);
	outBuf.append('\t');
	outBuf.append(_startPosStr);
	printOtherFields(outBuf);
}

void VcfRecord::print(QuickString &outBuf, const QuickString & start, const QuickString & end) const {
	outBuf.append(_chrName);
	outBuf.append('\t');
	outBuf.append(_startPosStr);
	printOtherFields(outBuf);

}

void VcfRecord::printNull(QuickString &outBuf) const {
	outBuf.append(".\t-1\t.\t.\t.\t-1");
	for (int i= startOtherIdx; i < _numPrintFields; i++) {
		outBuf.append("\t.");
	}
}

void VcfRecord::printOtherFields(QuickString &outBuf) const {
	outBuf.append('\t');
	outBuf.append(_name);
	outBuf.append('\t');
	outBuf.append(_varAlt);
	outBuf.append('\t');
	outBuf.append(_varRef);
	outBuf.append('\t');
	outBuf.append(_score);
	for (int i= 0; i < (int)_otherIdxs.size(); i++) {
		outBuf.append('\t');
		outBuf.append(*(_otherIdxs[i]));
	}

}

const QuickString &VcfRecord::getField(int fieldNum) const
{
	//a request for any of the first six fields will retrieve
	//chrom, start, end, name, score, and strand, in that order.
	//A request for field 6+ will go to the otherIdxs.

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
		return _varAlt;
		break;
	case 5:
		return _varRef;
		break;
	case 6:
		return _score;
	default:
		return BedPlusInterval::getField(fieldNum);
		break;
	}
}
