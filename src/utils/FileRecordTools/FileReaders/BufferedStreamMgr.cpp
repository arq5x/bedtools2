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
  	_eof(false)
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

	_mainBuf = new bufType[MAIN_BUF_READ_SIZE +1];
	memset(_mainBuf, 0, MAIN_BUF_READ_SIZE +1);

	return true;
}

bool BufferedStreamMgr::getTypeData()
{
	QuickString currScanBuffer;
	_inputStreamMgr->getSavedData(currScanBuffer);
	do {
		if (!_typeChecker.scanBuffer(currScanBuffer.c_str(), currScanBuffer.size()) && !_typeChecker.needsMoreData()) {
			return false;
		} else if (_typeChecker.needsMoreData()) {
			_inputStreamMgr->populateScanBuffer();
			currScanBuffer.clear();
			_inputStreamMgr->getSavedData(currScanBuffer);
		}
	} while (_typeChecker.needsMoreData());
	_inputStreamMgr->reset();
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
	memset(_mainBuf, 0, MAIN_BUF_READ_SIZE +1);

	_inputStreamMgr->getFinalStream()->read((char *)_mainBuf, MAIN_BUF_READ_SIZE);
	_mainBufCurrLen = _inputStreamMgr->getFinalStream()->gcount();
	_mainBufCurrStartPos = 0;
	return _mainBufCurrLen > 0;
}
