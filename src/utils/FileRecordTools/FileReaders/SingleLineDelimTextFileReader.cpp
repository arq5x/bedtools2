#include "SingleLineDelimTextFileReader.h"
#include <iostream>
#include <limits.h>
#include "BufferedStreamMgr.h"
#include "ParseTools.h"

bool abs_cmp(int i, int j) { return abs(i)<abs(j); }

SingleLineDelimTextFileReader::SingleLineDelimTextFileReader(int numFields, char delimChar)
: _numFields(numFields),
  _delimChar(delimChar),
  _fullHeaderFound(false),
  _currDataPos(0),
  _lineNum(0),
  _inheader(false)
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
		_lineNum++;
		if (!_bufStreamMgr->getLine(_sLine)) {
			return false;
		}
	}
	//after the first time we find a header, any other header line
	//doesn't count, and should not be added to the header, because it is a
	// "commented out record."
	_fullHeaderFound = true;

	//check to make sure line has something besides whitespace.
	bool hasNonSpace = false;
	int lineLen = (int)_sLine.size();
	for (int i=0; i < lineLen; i++) {
		if (!isspace(_sLine[i])) {
			hasNonSpace = true;
			break;
		}
	}
	if (!hasNonSpace) {
		return false;
	}

	if (wasHeader) {
		return true;
	}

	if (!findDelimiters()) {
		return false;
	}
	return true;
}


void SingleLineDelimTextFileReader::getField(int fieldNum, string &str) const {
	CHRPOS startPos = _delimPositions[fieldNum] +1;
	CHRPOS endPos = _delimPositions[fieldNum+1];
	str.assign(_sLine.c_str() + startPos, endPos - startPos);
}


void SingleLineDelimTextFileReader::getField(int fieldNum, CHRPOS &val) const {
	string temp;
	getField(fieldNum, temp);
	val = str2chrPos(temp.c_str());
}

void SingleLineDelimTextFileReader::getField(int fieldNum, char &val) const {
	val = _sLine[_delimPositions[fieldNum] +1];
}

void SingleLineDelimTextFileReader::appendField(int fieldNum, string &str) const {
	int startPos = _delimPositions[fieldNum] +1;
	int endPos = _delimPositions[fieldNum+1];
	str.append(_sLine.c_str() + startPos, endPos - startPos);
}

bool SingleLineDelimTextFileReader::detectAndHandleHeader()
{
	//not sure why the linker is giving me a hard time about
	//passing a non-const string to isHeaderLine, but
	//this const ref is a workaround.
	const string &sLine2 = _sLine;
	if (!isHeaderLine(sLine2) && (!(_inheader && _lineNum==1))) {
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
	int nfields = 0;
	const char* ptr = _sLine.c_str();
	const char* start = ptr;
	int target = _delimChar;

	while(ptr) {
		_delimPositions[nfields++] = ptr - start - 1;
		ptr = strchr(ptr, target);
		if(ptr) ptr ++;
	}
	_delimPositions[nfields] = _sLine.size();

	// count the number of fields allow for lines ending in \t\n
	if (nfields != _numFields) {
		cerr << "Error: line number " 
			 << _lineNum << " of file " 
			 << _filename << " has " 
			 << nfields << " fields, but " 
			 << _numFields 
			 << " were expected." 
			 << endl;
		exit(1);
	}
	return true;
}

CHRPOS SingleLineDelimTextFileReader::getVcfSVlen() {
	// tokenize the INFO field
	string info_str;
	vector<string> infofields;
	getField(VCF_TAG_FIELD, info_str);
	Tokenize(info_str, infofields, ';');

	// FOO=BAR;BIZ=BAM;SVLEN=100;END=200
	for (vector<string>::iterator f = infofields.begin(); f != infofields.end(); ++f) {
    	if (*f == ".") {
    		continue;
        }
        vector<string> keytoval;
        Tokenize(*f, keytoval, '=');
        //hey->val
        //SVLEN->100
        if (keytoval.size() == 2) {
        	if (keytoval.at(0) == "SVLEN") {
        		vector<CHRPOS> svlens;
        		Tokenize(keytoval.at(1), svlens, ',');
        		// are the multiple SVLENS?
        		if (svlens.size() == 1) {
        			return abs(svlens[0]);
        		}
        		else {
        			// return the abs_max SVLEN
        			CHRPOS max_len = *max_element(svlens.begin(),svlens.end(), abs_cmp);
        			return abs(max_len);
        		}
        	}
        	else if (keytoval.at(0) == "END") {
        		string start_str;
        		getField(1, start_str);
        		// length is END - POS + 1
        		return (int)(str2chrPos(keytoval.at(1)) - str2chrPos(start_str) + 1);
        	}
		}
    }
    // not found
    return INT_MIN;
}
