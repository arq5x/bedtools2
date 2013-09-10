/*
 * SingleLineTextFileReader.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef SINGLELINETEXTFILEREADER_H_
#define SINGLELINETEXTFILEREADER_H_

using namespace std;

#include "FileReader.h"
#include "QuickString.h"

class SingleLineDelimTextFileReader : public FileReader {
public:
	SingleLineDelimTextFileReader(int numFields, char delimChar = '\t');
	~SingleLineDelimTextFileReader();

	virtual bool readEntry();
	virtual int getNumFields() const { return _numFields; }
	virtual void getField(int numField, QuickString &val) const;
	virtual void getField(int numField, int &val); //this signaiture isn't const because it operates on an internal QuickString for speed.
	virtual void getField(int fieldNum, char &val) const;
	virtual void appendField(int fieldNum, QuickString &str) const;
	virtual const QuickString &getHeader() const { return _header; }
	virtual bool hasHeader() const { return _fullHeaderFound; }
protected:
	int _numFields;
	char _delimChar;
	QuickString _header;
	bool _fullHeaderFound;
	int _currDataPos;
	QuickString _sLine;
	int *_delimPositions;
	QuickString _currChromStr;
	QuickString _tempChrPosStr;
	int _lineNum;
	bool detectAndHandleHeader();
	bool findDelimiters();

};


#endif /* SINGLELINETEXTFILEREADER_H_ */
