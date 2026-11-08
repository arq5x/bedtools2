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
  _line(_sLine.c_str()),
  _lineLen(0),
  _lineNum(0),
  _inheader(false),
  _coordsOnly(false)
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
	if (!nextLine()) {
		return false;
	}
	_lineNum++;
	if (_lineLen == 0) {
		return false;
	}

	//scan the whole header in one call.
	bool wasHeader = false;
	while (detectAndHandleHeader()) { //header line
		_lineNum++;
		if (!nextLine()) {
			return false;
		}
	}
	//after the first time we find a header, any other header line
	//doesn't count, and should not be added to the header, because it is a
	// "commented out record."
	_fullHeaderFound = true;

	//check to make sure line has something besides whitespace.
	bool hasNonSpace = false;
	int lineLen = (int)_lineLen;
	for (int i=0; i < lineLen; i++) {
		if (!isspace(_line[i])) {
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


bool SingleLineDelimTextFileReader::nextLine() {
	return _bufStreamMgr->getLineRef(_sLine, _line, _lineLen);
}

void SingleLineDelimTextFileReader::getField(int fieldNum, string &str) const {
	CHRPOS startPos = _delimPositions[fieldNum] +1;
	CHRPOS endPos = _delimPositions[fieldNum+1];
	str.assign(_line + startPos, endPos - startPos);
}


void SingleLineDelimTextFileReader::getField(int fieldNum, CHRPOS &val) const {
	string temp;
	getPosField(fieldNum, val, temp);
}

bool SingleLineDelimTextFileReader::getPosField(int fieldNum, CHRPOS &val, string &str) const {
	const char *p = _line + _delimPositions[fieldNum] + 1;
	const char *e = _line + _delimPositions[fieldNum+1];
	const char *q = p;
	bool neg = (q < e && *q == '-');
	if (neg) q++;
	const char *digits = q;
	long long result = 0;
	while (q < e && (unsigned char)(*q - '0') <= 9) {
		result = result * 10 + (*q - '0');
		q++;
	}
	size_t nd = q - digits;
	// 1..18 digits cannot overflow; a leading zero or "-0" is not the plain form.
	if (q == e && nd >= 1 && nd <= 18 && (*digits != '0' || (nd == 1 && !neg))) {
		val = neg ? -result : result;
		return true;
	}
	str.assign(p, e - p);
	val = str2chrPos(str.c_str());
	return false;
}

void SingleLineDelimTextFileReader::getField(int fieldNum, char &val) const {
	val = _line[_delimPositions[fieldNum] +1];
}

void SingleLineDelimTextFileReader::appendField(int fieldNum, string &str) const {
	int startPos = _delimPositions[fieldNum] +1;
	int endPos = _delimPositions[fieldNum+1];
	str.append(_line + startPos, endPos - startPos);
}

// False only when isHeaderLine is certain to return false: its first character
// must be one of > ! # or the start of chrom, chr, browser, track or visibility,
// and "chr" followed by a digit is never one.
bool SingleLineDelimTextFileReader::mayBeHeader() const {
	unsigned char c = (unsigned char)_line[0];
	if (c == '>' || c == '!' || c == '#') {
		return true;
	}
	c |= 0x20;
	if (c == 'c') {
		return !(_lineLen >= 4 && (_line[1] | 0x20) == 'h' && (_line[2] | 0x20) == 'r' && isdigit((unsigned char)_line[3]));
	}
	return c == 'b' || c == 't' || c == 'v';
}

bool SingleLineDelimTextFileReader::detectAndHandleHeader()
{
	if (!(mayBeHeader() && isHeaderLine(string(_line, _lineLen))) && (!(_inheader && _lineNum==1))) {
		return false;
	}
	if (!_fullHeaderFound) {
		_header.append(_line, _lineLen);
		_header += "\n"; //add new line, since it was chomped by getline
	}
	return true;
}

bool SingleLineDelimTextFileReader::findDelimiters() {
	memset(_delimPositions, 0, (_numFields +1) * sizeof(int));
	//scan the line for delimiters, determine their positions
	int nfields = 0;
	const char* ptr = _line;
	const char* start = ptr;
	int target = _delimChar;

	while(ptr) {
		_delimPositions[nfields++] = ptr - start - 1;
		ptr = strchr(ptr, target);
		if(ptr) ptr ++;
	}
	_delimPositions[nfields] = (int)_lineLen;

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
