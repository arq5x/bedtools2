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
	// Same lines as getLine. When the line lies wholly inside the main buffer, ptr/len
	// point into it (NUL-terminated at ptr[len], valid until the next call); otherwise
	// the line is built in scratch. Returns false exactly when getLine would.
	bool getLineRef(string &scratch, const char *&ptr, size_t &len);
	// True when the line the last getLineRef returned lies in a chunk known to hold no NUL byte.
	bool lastLineNulFree() const { return _lineNulFree; }
	// Retain mode: every chunk read is kept (a fresh buffer per read) so that lines
	// handed out in place stay valid. Only for streams read in chunks of 8,191 bytes or more.
	void setRetain(bool val) { _retain = val && _useBufSize >= (size_t)GZIP_LINE_BUF_SIZE; }
	// True when the last getLineRef returned a NUL-free line in place in a retained chunk:
	// its chunk number (the position in the list releaseChunks returns) and byte offset.
	bool lastLineChunkPos(size_t &chunk, size_t &off) const {
		if (!_retain || !_lineNulFree) return false;
		chunk = _retained.size();
		off = _lineOff;
		return true;
	}
	// Hands every retained chunk, and the current one, to the caller, who now owns them.
	void releaseChunks(vector<char *> &out);
	BamTools::BamReader *getBamReader() { return _inputStreamMgr->getBamReader(); }
	static const int DEFAULT_MAIN_BUF_READ_SIZE = 1023;
	static const int PLAIN_FILE_BUF_READ_SIZE = 65535;
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
	bool _chunkNulFree;
	bool _lineNulFree;
	bool _retain;
	size_t _lineOff;
	vector<char *> _retained;
	string _currScanBuffer;

	//The minus ones in these constants are for leaving room for a null terminator after reading into buffers.
	static const int GZIP_LINE_BUF_SIZE = 8191; // 8K
	bool readFileChunk();
	bool getTypeData();
};


#endif /* BUFFEREDSTREAMMGR_H_ */
