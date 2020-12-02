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

bool BufferedStreamMgr::getLine(string &line)
{
	line.clear();

	if (_mainBufCurrStartPos >= _mainBufCurrLen) {
		if (!readFileChunk()) {
			_eof = true;
			return false;
		}
	}
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

	return retVal;
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
