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

BufferedStreamMgr::BufferedStreamMgr(const QuickString &filename)
: 	_inputStreamMgr(NULL),
	_mainBuf(NULL),
   	_filename(filename),
  	_mainBufCurrStartPos(0),
  	_mainBufCurrLen(0),
  	_eof(false),
  	_useBufSize(0),
  	_streamFinished(false)
{

}

BufferedStreamMgr::~BufferedStreamMgr()
{

	if (_mainBuf != NULL) {
		delete [] _mainBuf;
	}
	if (_inputStreamMgr != NULL) {
		delete _inputStreamMgr;
		_inputStreamMgr = NULL;
	}
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
		_typeChecker.setBam();
		return true;
	}
	if (!getTypeData()) {
		return false;
	}
	if (_inputStreamMgr->isGzipped()) {
		_useBufSize = GZIP_LINE_BUF_SIZE;
	} else {
		_useBufSize =  MAIN_BUF_READ_SIZE;
	}

	size_t trueBufSize = max(_useBufSize, (int)_currScanBuffer.size());
	_useBufSize = trueBufSize;
	_mainBuf = new bufType[_useBufSize +1];
	memset(_mainBuf, 0, _useBufSize +1);
	memcpy(_mainBuf, _currScanBuffer.c_str(), _currScanBuffer.size());
	_mainBufCurrLen = _currScanBuffer.size();

	return true;
}

bool BufferedStreamMgr::getTypeData()
{
	_currScanBuffer = _inputStreamMgr->getSavedData();
	_typeChecker.setFilename(_filename);
	do {
		if (!_typeChecker.scanBuffer(_currScanBuffer.c_str(), _currScanBuffer.size()) && !_typeChecker.needsMoreData()) {
			return false;
		} else if (_typeChecker.needsMoreData()) {
			_inputStreamMgr->populateScanBuffer();
			_currScanBuffer.append(_inputStreamMgr->getSavedData());
		}
	} while (_typeChecker.needsMoreData());
	if (_inputStreamMgr->resetStream()) {
		_currScanBuffer.clear();
	}
	return true;
}

bool BufferedStreamMgr::getLine(QuickString &line)
{
	line.clear();

	if (_mainBufCurrStartPos >= _mainBufCurrLen) {
		if (!readFileChunk()) {
			_eof = true;
			return false;
		}
	}
	while (1) {
		int searchPos = _mainBufCurrStartPos;
		while (searchPos < _mainBufCurrLen && _mainBuf[searchPos] != '\n') {
			searchPos++;
		}

		line.append((char *)_mainBuf + _mainBufCurrStartPos, searchPos - _mainBufCurrStartPos);
		_mainBufCurrStartPos = searchPos +1;
		if (searchPos == _mainBufCurrLen) { //hit end of buffer, but no newline yet
			if (!readFileChunk()) { //hit eof
				return true;
			}
		} else if (_mainBuf[searchPos] == '\n') {
			return true;
		}
	}
}

bool BufferedStreamMgr::readFileChunk()
{
	if (eof()) {
		return false;
	}
	memset(_mainBuf, 0, _useBufSize +1);
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
