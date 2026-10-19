/*
 * BufferedStreamMgr.cpp
 *
 *  Created on: Jul 9, 2013
 *      Author: nek3d
 */

#include "BufferedStreamMgr.h"
#include "CompressionTools.h"
#include "InputStreamMgr.h"
#include <fstream>
#include <stdint.h>

// getLine finds a line's end and its delimiters in one pass where SSE2 is
// available. Over the 64 bytes at s, bit k of d is set where s[k] == c, and
// bit k of e where s[k] is '\n' or '\0'.
#if defined(__SSE2__) && defined(__GNUC__)
#include <emmintrin.h>
#define LINE_SCAN_SSE2 1
static inline void scan64(const char *s, char c, uint64_t &d, uint64_t &e) {
	const __m128i cv = _mm_set1_epi8(c), nv = _mm_set1_epi8('\n'), zv = _mm_setzero_si128();
	d = e = 0;
	for (int k = 0; k < 64; k += 16) {
		__m128i v = _mm_loadu_si128((const __m128i *)(s + k));
		d |= (uint64_t)(unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(v, cv)) << k;
		e |= (uint64_t)(unsigned)_mm_movemask_epi8(_mm_or_si128(_mm_cmpeq_epi8(v, nv), _mm_cmpeq_epi8(v, zv))) << k;
	}
}
// Counts the delimiters whose bits are set in d, at offsets off + bit, in
// nf, storing delimiter j's offset to pos[j] while j <= cap.
static inline void storeDelims(uint64_t d, size_t off, int *pos, int cap, int &nf) {
	while (d) {
		if (nf <= cap) {
			pos[nf] = (int)(off + __builtin_ctzll(d));
		}
		nf++;
		d &= d - 1;
	}
}
#endif

BufferedStreamMgr::BufferedStreamMgr(const string &filename)
: 	_inputStreamMgr(NULL),
	_mainBuf(NULL),
   	_filename(filename),
  	_mainBufCurrStartPos(0),
  	_mainBufCurrLen(0),
  	_eof(false),
  	_useBufSize(DEFAULT_MAIN_BUF_READ_SIZE),
  	_streamFinished(false)
{

}

BufferedStreamMgr::~BufferedStreamMgr()
{

	delete [] _mainBuf;
	_mainBuf = NULL;

	delete _inputStreamMgr;
	_inputStreamMgr = NULL;
}

bool BufferedStreamMgr::init()
{
	_inputStreamMgr = new InputStreamMgr(_filename);
	if (!_inputStreamMgr->init()) {
		return false;
	}
	if (_inputStreamMgr->isBam()) {
		//there is a special check for a BAM file's magic string inside
		//the inputStreamMgr's init method. If it is found, we can
		//stop here.
		if(_inputStreamMgr->isCram())
			_typeChecker.setCram();
		else
			_typeChecker.setBam();
		return true;
	}
	if (!getTypeData()) {
		return false;
	}
	if (_inputStreamMgr->isGzipped()) {
		_useBufSize = GZIP_LINE_BUF_SIZE;
	}

	size_t trueBufSize = max(_useBufSize, _currScanBuffer.size());
	_useBufSize = (int)trueBufSize;
	_mainBuf = new bufType[_useBufSize +1];
	memset(_mainBuf, 0, _useBufSize +1);
	memcpy(_mainBuf, _currScanBuffer.c_str(), _currScanBuffer.size());
	_mainBufCurrLen = (int)_currScanBuffer.size();

	return true;
}

bool BufferedStreamMgr::getTypeData()
{
	_currScanBuffer = _inputStreamMgr->getSavedData();
	_typeChecker.setFilename(_filename);
	do {
		if (!_typeChecker.scanBuffer(_currScanBuffer.c_str(), _currScanBuffer.size(), _inputStreamMgr->getEofHit(), _inputStreamMgr->isCompressed()) && !_typeChecker.needsMoreData()) {
			return false;
		} else if (_typeChecker.needsMoreData()) {
			if (!_inputStreamMgr->populateScanBuffer()) {
				// We have found a file with just a header. If the file and record type are known,
				//break here. Otherwise, assign those types to empty file and record type, then break.
				if ((_typeChecker.getFileType() != FileRecordTypeChecker::UNKNOWN_FILE_TYPE) &&
					(_typeChecker.getRecordType() != FileRecordTypeChecker::UNKNOWN_RECORD_TYPE)) {
					break;
				} else {
					_typeChecker.setFileType(FileRecordTypeChecker::EMPTY_FILE_TYPE);
					_typeChecker.setRecordType(FileRecordTypeChecker::EMPTY_RECORD_TYPE);
					break;
				}
			}
			_currScanBuffer.append(_inputStreamMgr->getSavedData());
		}
	} while (_typeChecker.needsMoreData());
	if (_inputStreamMgr->resetStream()) {
		_currScanBuffer.clear();
	}
	return true;
}

bool BufferedStreamMgr::getLine(string &line, const char *&text, size_t &len,
	char delim, int *pos, int cap, int &nfields)
{
	nfields = 0;
	if (_mainBufCurrStartPos >= _mainBufCurrLen) {
		if (!readFileChunk()) {
			_eof = true;
			line.clear();
			text = line.c_str();
			len = 0;
			return false;
		}
	}
	// the line in place when its newline is in this chunk; no '\n' lies
	// before that one, so only '\r's are stripped
	char *start = (char *)_mainBuf + _mainBufCurrStartPos;
#ifdef LINE_SCAN_SSE2
	if (pos != NULL) {
		// 64 bytes a step while the step lies in the chunk. The bytes between
		// the line's end and its newline are '\r's, never delim, so the
		// delimiters before the newline are the line's own.
		size_t avail = (size_t)(_mainBufCurrLen - _mainBufCurrStartPos);
		int nf = 1;
		for (size_t off = 0; off + 64 <= avail; off += 64) {
			uint64_t dm, em;
			scan64(start + off, delim, dm, em);
			if (em) {
				size_t n = off + __builtin_ctzll(em);
				if (start[n] != '\n') {
					break;	// a NUL before the newline: the scan below
				}
				storeDelims(dm & (((uint64_t)1 << (n - off)) - 1), off, pos, cap, nf);
				_mainBufCurrStartPos += (int)n + 1;
				while (n > 0 && start[n - 1] == '\r') n--;
				start[n] = '\0';
				text = start;
				len = n;
				nfields = nf;
				return true;
			}
			storeDelims(dm, off, pos, cap, nf);
		}
	}
#endif
	char *eol = (char *)memchr(start, '\n', _mainBufCurrLen - _mainBufCurrStartPos);
	if (eol != NULL) {
		size_t n = (size_t)(eol - start);
		_mainBufCurrStartPos += (int)n + 1;
		while (n > 0 && start[n - 1] == '\r') n--;
		start[n] = '\0';
		text = start;
		len = n;
		return true;
	}

	line.clear();
	bool retVal = true;
	while (1) {
		int searchPos = _mainBufCurrLen;
		const char* eol_pos = (const char*)memchr((const char*)_mainBuf + _mainBufCurrStartPos, '\n', _mainBufCurrLen - _mainBufCurrStartPos);
		if(eol_pos) {
			searchPos = eol_pos - (const char*)_mainBuf;
		}

		line.append((char *)_mainBuf + _mainBufCurrStartPos, searchPos - _mainBufCurrStartPos);

		_mainBufCurrStartPos = searchPos +1;
		if (searchPos == _mainBufCurrLen) { //hit end of buffer, but no newline yet
			if (!readFileChunk()) { //hit eof
				retVal = true;
				break;
			}
		} else if (_mainBuf[searchPos] == '\n') {
			retVal = true;
			break;
		}
	}
	//strip any whitespace characters, such as DOS newline characters or extra tabs,
	//from the end of the line
	int lastPos = (int)line.size();
	while (lastPos > 0 && (line[lastPos-1] == '\n' || line[lastPos-1] == '\r')) lastPos--;
	line.resize(lastPos);

	text = line.c_str();
	len = line.size();
	return retVal;
}

size_t BufferedStreamMgr::readRaw(char *dst, size_t cap)
{
	size_t got = 0;
	if (_mainBufCurrStartPos < _mainBufCurrLen) {
		got = (size_t)(_mainBufCurrLen - _mainBufCurrStartPos);
		if (got > cap) {
			got = cap;
		}
		memcpy(dst, _mainBuf + _mainBufCurrStartPos, got);
		_mainBufCurrStartPos += (int)got;
		return got;
	}
	// readFileChunk's reads, without the chunk buffer: a short read ends
	// the stream, and so does an error
	if (_eof || _streamFinished || cap == 0) {
		_eof = _eof || cap != 0;
		return 0;
	}
	int n = _inputStreamMgr->read(dst, cap);
	if (n <= 0) {
		_streamFinished = true;
		_eof = true;
		return 0;
	}
	if ((size_t)n < cap) {
		_streamFinished = true;
	}
	return (size_t)n;
}

bool BufferedStreamMgr::readFileChunk()
{
	if (eof()) {
		return false;
	}
	_mainBufCurrStartPos = 0;

	if (!_streamFinished) {
		_mainBufCurrLen = _inputStreamMgr->read((char *)_mainBuf, _useBufSize);
		if (_mainBufCurrLen < _useBufSize) {
			_streamFinished = true;
		}
		return _mainBufCurrLen > 0;
	} else {
		_mainBufCurrLen = 0;
		return false;
	}
}
