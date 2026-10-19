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
	// The next line, without its newline and the '\r's before it, as text and
	// len, followed by a NUL. A line whose newline lies in the current chunk
	// stays where it is: the byte after it, a consumed '\n' or '\r', becomes
	// the NUL, and text is valid until the next call. Any other line (one
	// crossing a chunk boundary, or a last line without a newline) is copied
	// into line, and text points there.
	// Where SSE2 is available and pos is not NULL (delim is then none of
	// '\n', '\r' and '\0'), the pass that finds an in-place line's end also
	// finds its delimiters, as SingleLineDelimTextFileReader::findDelimiters
	// does: nfields gets the line's field count and the offset of delimiter j
	// (from 1) goes to pos[j] while j <= cap. Otherwise, and for a line that
	// holds a NUL, is copied into line, or ends past the chunk's last whole
	// 64 bytes from its start, nfields is 0 and pos may hold anything.
	bool getLine(string &line, const char *&text, size_t &len,
		char delim, int *pos, int cap, int &nfields);
	// The stream's bytes from where getLine would read next, for a reader
	// that splits lines itself: up to cap bytes into dst, what is left of
	// the current chunk first. Returns 0 at the end of the stream, where
	// getLine would return false, and eof() then holds.
	size_t readRaw(char *dst, size_t cap);
	// gzip, BGZF or BAM, whose reads can print errors as they go
	bool isCompressed() const { return _inputStreamMgr->isCompressed(); }
	BamTools::BamReader *getBamReader() { return _inputStreamMgr->getBamReader(); }
	static const int DEFAULT_MAIN_BUF_READ_SIZE = 65535; // 64K, minus one for a null terminator
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
