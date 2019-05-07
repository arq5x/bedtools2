/*
 * FileRecordTypeChecker.h
 *
 *  Created on: Nov 19, 2012
 *      Author: nek3d
 */

#ifndef FILERECORDTYPECHECKER_H_
#define FILERECORDTYPECHECKER_H_

#include <string>
#include <cstring>
#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <map>
#include <algorithm>
#include "PushBackStreamBuf.h"
#include "Tokenizer.h"

using namespace std;

class FileRecordTypeChecker {
public:

	FileRecordTypeChecker();
	~FileRecordTypeChecker() {}

	typedef enum  { UNKNOWN_FILE_TYPE, EMPTY_FILE_TYPE, SINGLE_LINE_DELIM_TEXT_FILE_TYPE,
			MULTI_LINE_ENTRY_TEXT_FILE_TYPE,
			GFF_FILE_TYPE, GZIP_FILE_TYPE, BAM_FILE_TYPE, VCF_FILE_TYPE} FILE_TYPE;

	typedef enum  { UNKNOWN_RECORD_TYPE, EMPTY_RECORD_TYPE, BED3_RECORD_TYPE, BED4_RECORD_TYPE, BEDGRAPH_RECORD_TYPE, BED5_RECORD_TYPE,
		BED6_RECORD_TYPE, BED12_RECORD_TYPE, BED_PLUS_RECORD_TYPE, BED6_PLUS_RECORD_TYPE, BAM_RECORD_TYPE, VCF_RECORD_TYPE, GFF_RECORD_TYPE,
		GFF_PLUS_RECORD_TYPE, NO_POS_PLUS_RECORD_TYPE} RECORD_TYPE;

	void setFilename(const string & filename) { _filename = filename; }
	bool scanBuffer(const char *buf, size_t len, bool eofHit, bool isCompressed = false);
	bool needsMoreData() const { return _insufficientData; }

	bool recordTypeHasName(RECORD_TYPE type) const { return _hasName.find(type) != _hasName.end(); }
	bool recordTypeHasScore(RECORD_TYPE type) const { return _hasScore.find(type) != _hasScore.end(); }
	bool recordTypeHasStrand(RECORD_TYPE type) const { return _hasStrand.find(type) != _hasStrand.end(); }

	FILE_TYPE getFileType() const { return _fileType; }
	void setFileType(FILE_TYPE type) { _fileType = type; }


	RECORD_TYPE getRecordType() const { return _recordType; }
	void setRecordType(RECORD_TYPE type) { _recordType = type; }

	const string &getFileTypeName() const {
		return _fileTypeNames.find(_fileType)->second;
	}
	const string &getRecordTypeName() const {
		return _recordTypeNames.find(_recordType)->second;
	}

	bool isBinary() const { return _isBinary; }
	bool isBam() const { return _isBAM; }
	bool isCram() const { return _isCRAM; }
	bool isGzipped() const { return _isGzipped; }

	void setBam(); //call only if you're SURE the file is BAM!
	void setCram();
	void setIsGroupBy(bool val) { _isGroupBy = val; } // When using groupBy,


	bool isText() const { return _isText; }
	bool isDelimited() const { return _isDelimited; }
	char getDelimChar() const { return _delimChar; }
	int getNumFields() const { return _numFields; }

	bool isVcf() const;
	bool isBed() const { return _isBed; }
	bool isBed3() const { return (_isBed && _numFields == 3); }
	bool isBed4() const { return (_isBed && _numFields == 4 && !_fourthFieldNumeric); }
	bool isBedGraph() const { return (_isBed && _numFields == 4 && _fourthFieldNumeric); }
	bool isBed5() const { return (_isBed && _numFields == 3); }
	bool isBed6() const { return (_isBed && _numFields == 6); }
	bool isBedPlus() const { return (_isBed && _numFields > 6 && _numFields != 12); }
	bool isBed12() const { return (_isBed && _numFields == 12); }
	bool isGFF() const { return _isGFF; }

	void setInHeader(bool val) { _inheader = val; }




private:
	FILE_TYPE _fileType;
	RECORD_TYPE _recordType;

	string _filename; //useful for reporting errors with file.
	Tokenizer _tokenizer;

	int _firstValidDataLineIdx;
	int _numBytesInBuffer; //this will hold the length of the buffer after the scan.

	int _numFields;
	bool _isBinary;
	bool _isText;
	bool _isBed;
	bool _isDelimited;
	char _delimChar;
	bool _isVCF;
	bool _isBAM;
	bool _isCRAM;
	bool _isGFF;
	bool _isGFFplus;
	bool _isGzipped;
	bool _isCompressed;
	bool _insufficientData; //set to true if scan buffer had only header lines.
	bool _fourthFieldNumeric; //this is just to distinguish between Bed4 and BedGraph files.
	bool _givenEmptyBuffer;
	bool _isGroupBy;

	map<RECORD_TYPE, string> _recordTypeNames;
	map<FILE_TYPE, string> _fileTypeNames;

	map<RECORD_TYPE, bool> _hasName;
	map<RECORD_TYPE, bool> _hasScore;
	map<RECORD_TYPE, bool> _hasStrand;

	bool _eofHit;
	bool _inheader;

	bool isBinaryBuffer(const char *buffer, size_t len);
	bool isBAM(const char *buffer);
	bool handleTextFormat(const char *buffer, size_t len);
	bool isTextDelimtedFormat(const char *buffer, size_t len);
	bool isBedFormat();
	bool isVCFformat(const char *buffer);
	bool isGFFformat();
	bool delimiterTesting(vector<int> &counts, char suspectChar);
	bool isStrandField(int field);
	bool passesBed5();
	bool passesBed6();
	bool passesBed12();

};

#endif /* FILERECORDTYPECHECKER_H_ */
