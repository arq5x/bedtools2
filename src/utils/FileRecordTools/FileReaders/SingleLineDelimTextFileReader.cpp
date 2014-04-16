#include "SingleLineDelimTextFileReader.h"
#include <iostream>
#include "BufferedStreamMgr.h"
#include "ParseTools.h"

SingleLineDelimTextFileReader::SingleLineDelimTextFileReader(int numFields, char delimChar)
: _numFields(numFields),
  _delimChar(delimChar),
  _fullHeaderFound(false),
  _currDataPos(0),
  _lineNum(0)
{
	_delimPositions = new int[numFields +1];
}

SingleLineDelimTextFileReader::~SingleLineDelimTextFileReader()
{
	delete [] _delimPositions;
	_delimPositions = NULL;
}

bool SingleLineDelimTextFileReader::readEntry()
{
	if (!_isFileOpen) {
		return false;
	}

	if (_bufStreamMgr->eof()) {
		return false;
	}
	if (!_bufStreamMgr->getLine(_sLine)) {
		return false;
	}
	_lineNum++;
	if (_sLine.empty()) {
		return false;
	}

	//scan the whole header in one call.
	bool wasHeader = false;
	while (detectAndHandleHeader()) { //header line
		if (!_bufStreamMgr->getLine(_sLine)) {
			return false;
			_lineNum++;
		}
	}
	//after the first time we find a header, any other header line
	//doesn't count, and should not be added to the header, because it is a
	// "commented out record."
	_fullHeaderFound = true;

	//check to make sure line has something besides whitespace.
	bool hasNonSpace = false;
	int lineLen = _sLine.size();
	for (int i=0; i < lineLen; i++) {
		if (!isspace(_sLine[i])) {
			hasNonSpace = true;
			break;
		}
	}
	if (!hasNonSpace) {
		return false;
	}

	//trim off any white space from end of line.
	int currPos = lineLen-1;
	while (isspace(_sLine[currPos])) {
		currPos--;
	}
	_sLine.resize(currPos +1);

	if (wasHeader) {
		return true;
	}

	if (!findDelimiters()) {
		return false;
	}
	return true;
}


void SingleLineDelimTextFileReader::getField(int fieldNum, QuickString &str) const {
	int startPos = _delimPositions[fieldNum] +1;
	int endPos = _delimPositions[fieldNum+1];
	str.assign(_sLine.c_str() + startPos, endPos - startPos);
}


void SingleLineDelimTextFileReader::getField(int fieldNum, int &val) {
	getField(fieldNum, _tempChrPosStr);
	val = str2chrPos(_tempChrPosStr.c_str());
}

void SingleLineDelimTextFileReader::getField(int fieldNum, char &val) const {
	val = _sLine[_delimPositions[fieldNum] +1];
}

void SingleLineDelimTextFileReader::appendField(int fieldNum, QuickString &str) const {
	int startPos = _delimPositions[fieldNum] +1;
	int endPos = _delimPositions[fieldNum+1];
	str.append(_sLine.c_str() + startPos, endPos - startPos);
}

bool SingleLineDelimTextFileReader::detectAndHandleHeader()
{
	//not sure why the linker is giving me a hard time about
	//passing a non-const QuickString to isHeaderLine, but
	//this const ref is a workaround.
	const QuickString &sLine2 = _sLine;
	if (!isHeaderLine(sLine2)) {
		return false;
	}
	if (!_fullHeaderFound) {
		_header += _sLine;
		_header += "\n"; //add new line, since it was chomped by getline
	}
	return true;
}

bool SingleLineDelimTextFileReader::findDelimiters() {
	memset(_delimPositions, 0, (_numFields +1) * sizeof(int));
	//scan the line for delimiters, determine their positions
	_delimPositions[0] = -1;
	int currField=1;
	int len = (int)_sLine.size();
	for (int i=0; i < len; i++) {
		if (_sLine[i] == _delimChar) {
			_delimPositions[currField] = i;
			currField++;
		}
	}
	_delimPositions[currField] = len;
	if (currField != _numFields) {
		cerr << "Error: line number " << _lineNum << " of file " << _filename << " has " << currField << " fields, but " << _numFields << " were expected." << endl;
		exit(1);
	}
	return true;
}
