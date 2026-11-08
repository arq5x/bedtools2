/*
 * SingleLineTextFileReader.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef SINGLELINETEXTFILEREADER_H_
#define SINGLELINETEXTFILEREADER_H_

#include <algorithm>
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
	int readCoordsLine(const string &lastChr, CHRPOS &s, CHRPOS &e);
	virtual int getNumFields() const { return _numFields; }
	const char *getLineData() const { return _line; }
	size_t getLineLength() const { return _lineLen; }
	const int *getDelimPositions() const { return _delimPositions; }
	// True and val set only for a non-negative plain decimal (no leading zero, 1-18 digits).
	bool getPlainPos(int fieldNum, CHRPOS &val) const;
	virtual void getField(int numField, string &val) const;
	virtual void getField(int numField, CHRPOS &val) const; //this signaiture isn't const because it operates on an internal string for speed.
	virtual void getField(int fieldNum, char &val) const;
	// Parses the field as str2chrPos would, without copying it. Returns true when the
	// field's text is exactly the plain decimal form of val (str is then untouched);
	// otherwise str receives the field text.
	bool getPosField(int fieldNum, CHRPOS &val, string &str) const;
	// When set, Bed6/BedPlus records parse only chrom, start and end.
	void setCoordsOnly(bool val) { _coordsOnly = val; }
	bool getCoordsOnly() const { return _coordsOnly; }
	virtual void appendField(int fieldNum, string &str) const;
	virtual const string &getHeader() const { return _header; }
	virtual bool hasHeader() const { return _fullHeaderFound; }
	virtual void setInHeader(bool val) { _inheader = val; }

protected:
	int _numFields;
	char _delimChar;
	string _header;
	bool _fullHeaderFound;
	int _currDataPos;
	string _sLine; // scratch for lines not wholly inside the stream buffer
	const char *_line;
	size_t _lineLen;
	bool nextLine();
	bool finishEntry();
	int *_delimPositions;
	string _currChromStr;
	string _tempChrPosStr;
	int _lineNum;
	bool _inheader;
	bool _coordsOnly;
	bool detectAndHandleHeader();
	bool mayBeHeader() const;
	bool findDelimiters();

	//This is actually a very specialized function strictly for VCF
	//records to read and process extra information about Structural Variants.
	static const int VCF_TAG_FIELD = 7;
	CHRPOS getVcfSVlen();


};


#endif /* SINGLELINETEXTFILEREADER_H_ */
