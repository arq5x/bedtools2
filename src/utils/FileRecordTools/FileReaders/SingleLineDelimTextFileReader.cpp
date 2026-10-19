#include "SingleLineDelimTextFileReader.h"
#include <iostream>
#include <limits.h>
#include <stdint.h>
#include "BufferedStreamMgr.h"
#include "ParseTools.h"

// findDelimiters reads a line 16 bytes at a time where SSE2 is available
// and one byte at a time elsewhere; both stop at the first NUL, as strchr does.
#if defined(__SSE2__) && defined(__GNUC__)
#include <emmintrin.h>
#define DELIM_SCAN_SSE2 1
// Bit k of d is set where s[k] == c, and bit k of z where s[k] == '\0', k < 16.
static inline void scan16(const char *s, char c, unsigned &d, unsigned &z) {
	__m128i v = _mm_loadu_si128((const __m128i *)s);
	d = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, _mm_set1_epi8(c)));
	z = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, _mm_setzero_si128()));
}
#endif

// The number of fields in the NUL-terminated line s split at delim. The
// offset of delimiter j (from 1) goes to pos[j] while j <= cap; delimiters
// past cap are counted but not stored.
static int scanDelimsBytes(const char *s, char delim, int *pos, int cap) {
	int nfields = 1;
	for (int i = 0; ; i++) {
		char c = s[i];
		if (c == delim) {
			if (nfields <= cap) {
				pos[nfields] = i;
			}
			nfields++;
		} else if (c == '\0') {
			return nfields;
		}
	}
}

#ifdef DELIM_SCAN_SSE2
// As scanDelimsBytes, for a line of n >= 16 bytes counting its terminator.
// Each block of up to 64 bytes becomes one mask of its delimiters and one of
// its NULs; a block's last partial 16 bytes are read by a load ending at the
// terminator, which overlaps bytes already read.
static int scanDelimsSSE2(const char *s, size_t n, char delim, int *pos, int cap) {
	int nfields = 1;
	for (size_t i = 0; ; i += 64) {
		size_t blk = n - i < 64 ? n - i : 64;
		size_t full = blk & ~(size_t)15;
		uint64_t dm = 0, zm = 0;
		unsigned d, z;
		for (size_t k = 0; k < full; k += 16) {
			scan16(s + i + k, delim, d, z);
			dm |= (uint64_t)d << k;
			zm |= (uint64_t)z << k;
		}
		size_t r = blk - full;
		if (r) {
			scan16(s + n - 16, delim, d, z);
			dm |= (uint64_t)(d >> (16 - r)) << full;
			zm |= (uint64_t)(z >> (16 - r)) << full;
		}
		if (zm) {
			dm &= (zm & (0 - zm)) - 1;	// only the delimiters before the first NUL
		}
		while (dm) {
			if (nfields <= cap) {
				pos[nfields] = (int)(i + __builtin_ctzll(dm));
			}
			nfields++;
			dm &= dm - 1;
		}
		if (zm) {
			return nfields;
		}
	}
}
#endif

bool abs_cmp(int i, int j) { return abs(i)<abs(j); }

SingleLineDelimTextFileReader::SingleLineDelimTextFileReader(int numFields, char delimChar)
: _numFields(numFields),
  _delimChar(delimChar),
  _fullHeaderFound(false),
  _currDataPos(0),
  _lineNum(0),
  _inheader(false),
  _takeAfterHeader(false)
{
	_delimPositions = new int[numFields +1];
	_scanPositions = new int[numFields +1];
	_scanFields = 0;
	_scanInGetLine = delimChar != '\n' && delimChar != '\r' && delimChar != '\0';
	_line = _sLine.c_str();
	_lineLen = 0;
}

SingleLineDelimTextFileReader::~SingleLineDelimTextFileReader()
{
	delete [] _delimPositions;
	_delimPositions = NULL;
	delete [] _scanPositions;
	_scanPositions = NULL;
}

bool SingleLineDelimTextFileReader::readEntry()
{
	if (!_isFileOpen) {
		return false;
	}

	if (_bufStreamMgr->eof()) {
		return false;
	}
	if (!_bufStreamMgr->getLine(_sLine, _line, _lineLen, _delimChar,
			_scanInGetLine ? _scanPositions : NULL, _numFields, _scanFields)) {
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
		if (!_bufStreamMgr->getLine(_sLine, _line, _lineLen, _delimChar,
				_scanInGetLine ? _scanPositions : NULL, _numFields, _scanFields)) {
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


// readEntry, one line at a time: readEntry's first line returns early when
// empty, while a line its header loop reads goes straight to the header test.
bool SingleLineDelimTextFileReader::takeLine(const char *line, size_t len)
{
	_lineNum++;
	if (!_takeAfterHeader && len == 0) {
		return false;
	}
	_line = line;
	_lineLen = len;
	_scanFields = 0;
	if (detectAndHandleHeader()) {
		_takeAfterHeader = true;
		return false;
	}
	_takeAfterHeader = false;
	_fullHeaderFound = true;

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
	return findDelimiters();
}

void SingleLineDelimTextFileReader::setLine(const char *line, size_t len) {
	_sLine.assign(line, len);
	_line = _sLine.c_str();
	_lineLen = _sLine.size();
	_scanFields = 0;
	findDelimiters();
}

void SingleLineDelimTextFileReader::getField(int fieldNum, string &str) const {
	CHRPOS startPos = _delimPositions[fieldNum] +1;
	CHRPOS endPos = _delimPositions[fieldNum+1];
	str.assign(_line + startPos, endPos - startPos);
}


void SingleLineDelimTextFileReader::getField(int fieldNum, CHRPOS &val) const {
	string temp;
	getField(fieldNum, temp);
	val = str2chrPos(temp.c_str());
}

void SingleLineDelimTextFileReader::getField(int fieldNum, char &val) const {
	val = _line[_delimPositions[fieldNum] +1];
}

void SingleLineDelimTextFileReader::appendField(int fieldNum, string &str) const {
	int startPos = _delimPositions[fieldNum] +1;
	int endPos = _delimPositions[fieldNum+1];
	str.append(_line + startPos, endPos - startPos);
}

bool SingleLineDelimTextFileReader::detectAndHandleHeader()
{
	if (!isHeaderLine(_line, _lineLen) && (!(_inheader && _lineNum==1))) {
		return false;
	}
	if (!_fullHeaderFound) {
		_header.append(_line, _lineLen);
		_header += "\n"; //add new line, since it was chomped by getline
	}
	return true;
}

bool SingleLineDelimTextFileReader::findDelimiters() {
	//scan the line for delimiters, determine their positions: field i runs
	//from _delimPositions[i] + 1 to _delimPositions[i + 1]. Only the
	//_numFields + 1 entries are written; a line with more fields is counted.
	const char *line = _line;
	size_t len = _lineLen;
	int nfields;
	if (_scanFields > 0) {
		// getLine found them
		nfields = _scanFields;
		_scanFields = 0;
		int *found = _scanPositions;
		_scanPositions = _delimPositions;
		_delimPositions = found;
		_delimPositions[0] = -1;
	} else {
		_delimPositions[0] = -1;
#ifdef DELIM_SCAN_SSE2
		if (len >= 15) {
			nfields = scanDelimsSSE2(line, len + 1, _delimChar, _delimPositions, _numFields);
		} else {
			nfields = scanDelimsBytes(line, _delimChar, _delimPositions, _numFields);
		}
#else
		nfields = scanDelimsBytes(line, _delimChar, _delimPositions, _numFields);
#endif
	}

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
	_delimPositions[nfields] = len;
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
