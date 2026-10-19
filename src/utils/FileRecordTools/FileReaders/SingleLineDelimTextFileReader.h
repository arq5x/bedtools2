/*
 * SingleLineTextFileReader.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef SINGLELINETEXTFILEREADER_H_
#define SINGLELINETEXTFILEREADER_H_

#include <algorithm>
#include <cstring>
#include "FileReader.h"
#include "string.h"
#include "lineFileUtilities.h"

typedef int64_t CHRPOS;

class SingleLineDelimTextFileReader : public FileReader {
public:
	//Allow VCF records to access a specialized private method.
	//See end of class declaration for details.
	friend class VcfRecord;

	SingleLineDelimTextFileReader(int numFields, char delimChar = '\t');
	~SingleLineDelimTextFileReader();

	virtual bool readEntry();
	virtual int getNumFields() const { return _numFields; }
	char getDelimChar() const { return _delimChar; }
	virtual void getField(int numField, string &val) const;
	virtual void getField(int numField, CHRPOS &val) const; //this signaiture isn't const because it operates on an internal string for speed.
	virtual void getField(int fieldNum, char &val) const;
	virtual void appendField(int fieldNum, string &str) const;
	// Field fieldNum's bytes in the current line, as getField copies them.
	const char *fieldStart(int fieldNum) const { return _line + _delimPositions[fieldNum] + 1; }
	size_t fieldLength(int fieldNum) const { return (size_t)(_delimPositions[fieldNum + 1] - _delimPositions[fieldNum] - 1); }
	virtual const string &getHeader() const { return _header; }
	virtual bool hasHeader() const { return _fullHeaderFound; }
	virtual void setInHeader(bool val) { _inheader = val; }

	// The line readEntry last read (its bytes, followed by a NUL, valid
	// until the next read), and a way to parse a copy of one held
	// elsewhere: setLine finds its delimiters as readEntry did, so getField
	// then returns the same fields. The line must have the field count
	// readEntry checked, since setLine reports no error of its own.
	const char *currentLine() const { return _line; }
	size_t currentLineLength() const { return _lineLen; }
	void setLine(const char *line, size_t len);
	// Whether the line readEntry last read holds a NUL. findDelimiters
	// found all its fields before the first NUL, so one can lie only in
	// the last field.
	bool lineHasNul() const {
		size_t from = (size_t)(_delimPositions[_numFields - 1] + 1);
		return memchr(_line + from, '\0', _lineLen - from) != NULL;
	}

	// Reading lines split from the stream elsewhere (BinTree's parallel
	// load), from the start of the file: each line is given in order, as
	// getLine would have returned it (its bytes up to the newline, '\r's
	// before it stripped, a NUL after them). takeLine does what readEntry
	// does with such a line, with the same header handling and errors, and
	// returns whether it is an entry; the line stays current, in place,
	// until the next. skipLines(n) stands for n lines that takeLine would
	// have read as entries, none of them a header or blank.
	bool takeLine(const char *line, size_t len);
	void skipLines(size_t n) {
		_lineNum += (int)n;
		_fullHeaderFound = true;
		_takeAfterHeader = false;
	}
	// The current line's delimiter positions (getNumFields() + 1 of them),
	// and a way to make current a line kept elsewhere with the positions
	// found for it then, as readEntry left them (IntersectFile's batches).
	const int *delimPositions() const { return _delimPositions; }
	void showLine(const char *line, size_t len, const int *delimPositions) {
		_line = line;
		_lineLen = len;
		memcpy(_delimPositions, delimPositions, (size_t)(_numFields + 1) * sizeof(int));
	}
	// how many lines readEntry and takeLine have read, and whether the
	// first line is a header whatever it holds (-header)
	int linesRead() const { return _lineNum; }
	bool firstLineIsHeader() const { return _inheader; }

protected:
	int _numFields;
	char _delimChar;
	string _header;
	bool _fullHeaderFound;
	int _currDataPos;
	// The current line, _lineLen bytes and a NUL at _line: in the stream's
	// buffer where getLine left it there, else in _sLine.
	const char *_line;
	size_t _lineLen;
	string _sLine;
	int *_delimPositions;
	// Where getLine found the current line's delimiters in the pass that
	// found its end (_scanFields > 0, their count), findDelimiters swaps
	// _scanPositions, where getLine stored them, in for _delimPositions.
	int *_scanPositions;
	int _scanFields;
	bool _scanInGetLine;
	string _currChromStr;
	string _tempChrPosStr;
	int _lineNum;
	bool _inheader;
	// takeLine's last line was a header, so that readEntry's header loop
	// would be reading the next one
	bool _takeAfterHeader;
	bool detectAndHandleHeader();
	bool findDelimiters();

	//This is actually a very specialized function strictly for VCF
	//records to read and process extra information about Structural Variants.
	static const int VCF_TAG_FIELD = 7;
	CHRPOS getVcfSVlen();


};


#endif /* SINGLELINETEXTFILEREADER_H_ */
