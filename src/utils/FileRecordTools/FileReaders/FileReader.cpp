#include <iostream>
#include <cstring>
#include "FileReader.h"
#include "BufferedStreamMgr.h"

FileReader::FileReader()
:
//  _inputStream(NULL),
  _bufStreamMgr(NULL),
  _isFileOpen(false),
//  _mustDeleteInputStream(false),
//  _externalInputStream(false),
//  _useBufStream(true),
  _currChromId(-1)
{
}

FileReader::~FileReader() {
}

bool FileReader::open() {

	if (_isFileOpen) {
		return true;
	} else {
		fprintf(stderr, "Error: bad input stream.\n");
		exit(1);
	}
	return false; //can't reach this point, but it's here to satisfy the compiler.
//
//	printf("Inside FileReader::open.\n");
//	if (!_externalInputStream && _inputStream == NULL) {
//		_inputStream = new ifstream(_filename.c_str(), ios::in);
//		_mustDeleteInputStream = true;
//	}
//	if (_inputStream == NULL || !_inputStream->good()) {
//		fprintf(stderr, "Error: bad input stream.\n");
//		exit(1);
//	}
//
//	_isFileOpen = true;
//	return true;
}

void FileReader::close() {
//	if (_mustDeleteInputStream) {
//		delete _inputStream;
//	}
//	return;
}

bool FileReader::eof() const {
//	return _useBufStream ? _bufStreamMgr->eof() : _inputStream->eof();
	return _bufStreamMgr->eof();
}
