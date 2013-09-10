/*
 * BufferedStreamMgr.h
 *
 *  Created on: Jul 9, 2013
 *      Author: nek3d
 */

#ifndef BUFFEREDSTREAMMGR_H_
#define BUFFEREDSTREAMMGR_H_

using namespace std;

#include <iostream>
#include "QuickString.h"
#include "FileRecordTypeChecker.h"
#include "InputStreamMgr.h"

class BufferedStreamMgr {
public:
	BufferedStreamMgr(const QuickString &filename);
	~BufferedStreamMgr();

	bool init();

	const FileRecordTypeChecker & getTypeChecker() const { return _typeChecker; }
	istream *getStream() { return _inputStreamMgr->getFinalStream(); }

	bool eof() const { return _eof; }
	bool getLine(QuickString &line);

private:
	InputStreamMgr *_inputStreamMgr;
	typedef unsigned char bufType;
	bufType *_mainBuf;

	FileRecordTypeChecker _typeChecker;
	QuickString _filename;

	int _mainBufCurrStartPos;
	int _mainBufCurrLen;
	bool _eof;
	//The minus ones in these constants are for leaving room for a null terminator after reading into buffers.
	static const int MAIN_BUF_READ_SIZE = 67108863; //64 Mb minus 1
	static const int TYPE_CHECK_READ_SIZE = 4095; // 4K
	static const int GZIP_LINE_BUF_SIZE = 16384; //16K
	bool readFileChunk();
	bool getTypeData();
//	void resetStream();
};


#endif /* BUFFEREDSTREAMMGR_H_ */
