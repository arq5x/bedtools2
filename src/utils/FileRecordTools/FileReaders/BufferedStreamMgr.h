/*
 * BufferedStreamMgr.h
 *
 *  Created on: Jul 9, 2013
 *      Author: nek3d
 */

#ifndef BUFFEREDSTREAMMGR_H_
#define BUFFEREDSTREAMMGR_H_

#include <iostream>
#include "string.h"
#include "FileRecordTypeChecker.h"
#include "InputStreamMgr.h"

class BufferedStreamMgr {
public:
	BufferedStreamMgr(const string &filename);
	~BufferedStreamMgr();

	bool init();

	FileRecordTypeChecker & getTypeChecker() { return _typeChecker; }

	bool eof() const { return _eof; }
	bool getLine(string &line);
	BamTools::BamReader *getBamReader() { return _inputStreamMgr->getBamReader(); }
	static const int DEFAULT_MAIN_BUF_READ_SIZE = 1023;
	void setIoBufSize(size_t val) { _useBufSize = val; }
private:
	InputStreamMgr *_inputStreamMgr;
	typedef unsigned char bufType;
	bufType *_mainBuf;

	FileRecordTypeChecker _typeChecker;
	string _filename;

	int _mainBufCurrStartPos;
	int _mainBufCurrLen;
	bool _eof;
	size_t _useBufSize;
	bool _streamFinished;
	string _currScanBuffer;

	//The minus ones in these constants are for leaving room for a null terminator after reading into buffers.
	static const int GZIP_LINE_BUF_SIZE = 8191; // 8K
	bool readFileChunk();
	bool getTypeData();
};


#endif /* BUFFEREDSTREAMMGR_H_ */
