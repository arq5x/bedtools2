/*
 * InputStreamMgr.cpp
 *
 *  Created on: Mar 21, 2013
 *      Author: nek3d
 */

#include "InputStreamMgr.h"
#include <cstring> //for memset
#include "gzstream.h"
#include "CompressionTools.h"

InputStreamMgr::InputStreamMgr(const QuickString &filename, bool buildScanBuffer)
:
 _filename(filename),
 _pushBackStreamBuf(NULL),
 _inputFileStream(NULL),
 _infStreamBuf(NULL),
 _oldInputStream(NULL),
 _isStdin(false),
 _isGzipped(false),
 _isBam(false),
 _numBytesInBuffer(0)
{
	_possibleBamCode.resize(4, 0);
}


InputStreamMgr::~InputStreamMgr() {
	if (_pushBackStreamBuf != NULL) {
		delete _pushBackStreamBuf;
		_pushBackStreamBuf = NULL;
	}
	if (_inputFileStream != NULL) {
		delete _inputFileStream;
		_inputFileStream = NULL;
	}
	if (_oldInputStream != NULL) {
		delete _oldInputStream;
		_oldInputStream = NULL;
	}
	if (_infStreamBuf != NULL) {
		delete _infStreamBuf;
		_infStreamBuf = NULL;
	}
}

bool InputStreamMgr::init()
{
	if (_filename == "-" || _filename == "stdin") { //stdin
		_isStdin = true;
		//peek at the first char of stdin to see if this is gzipped.
		if ((unsigned char)cin.peek()  == 0x1f) {
			_isGzipped = true;
		}
		_pushBackStreamBuf = new PushBackStreamBuf(cin.rdbuf());
	} else {
		_inputFileStream = new ifstream(_filename.c_str());
		if (_inputFileStream->fail()) {
			cerr << "Error: Unable to open file " << _filename << ". Exiting." << endl;
			delete _inputFileStream;
			_inputFileStream = NULL;
			exit(1);
		}
		//peek at the first char of stdin to see if this is gzipped.
		if ((unsigned char)_inputFileStream->peek()  == 0x1f) {
			_isGzipped = true;
		}
		_pushBackStreamBuf = new PushBackStreamBuf(_inputFileStream->rdbuf());
	}
	//now we have a PushBackStreamBuf. Make a new stream.
	_finalInputStream = new istream(_pushBackStreamBuf);
	populateScanBuffer();
	return true;
}

void InputStreamMgr::populateScanBuffer()
{
	_scanBuffer.clear();
	int numChars=0;
	int currChar = 0;
	while (1) {
		currChar = _pushBackStreamBuf->sbumpc();
		//Stop when EOF hit.
		if (currChar == EOF) {
			break;
		}
		numChars++;
		_scanBuffer.push_back(currChar);
		if (_isGzipped) {
			if (bamDetected(numChars, currChar)) {
				return;
			}
			_compressedSaveData.push_back(currChar);
		}

		//For non-gzip, also stop if we have the minimum number of bytes and newline is hit.
		//For gzip, stop at SCAN_BUFFER_SIZE.
		if ((!_isGzipped && (currChar == '\n' && numChars >= MIN_SCAN_BUFFER_SIZE )) || (_isGzipped && numChars >= SCAN_BUFFER_SIZE)) {
			break;
		}
	}
	_numBytesInBuffer = _scanBuffer.size();

	//append it to the savedDataStr. If it's gzipped, decompress it first.
	if (_isGzipped) {
		decompressBuffer();
	} else {
		_scanBuffer.toStr(_saveDataStr, true);
	}
}

bool InputStreamMgr::bamDetected(int numChars, int currChar)
{
	//Look for the BAM magic string "BAM\1" in the first fouur characters of the input stream.
	//In compressed form, the first char is the gzip signifier, which was already found.
	//The next three are the integers 139, 8, and 4.
	if (numChars < 5) {
		_possibleBamCode[numChars -1] = currChar;
		//special: test for BAM
		if (numChars == 4 && _possibleBamCode[1] == 139 && _possibleBamCode[2] == 8 && _possibleBamCode[3] == 4) {
			//BAM detected.
			_pushBackStreamBuf->pushBack(_scanBuffer);
			_isBam = true;
			_numBytesInBuffer = 4;
			return true;
		}
	}
	return false;
}

void InputStreamMgr::decompressBuffer()
{
	//allocate an array to hold uncompressed data.
	uInt maxDecompressSize = 20 * _numBytesInBuffer;
	unsigned char *newScanBuffer = new unsigned char[maxDecompressSize];
	memset(newScanBuffer, 0, maxDecompressSize);

	unsigned int numDecompressChars = inflateGzippedArray(_scanBuffer, newScanBuffer, maxDecompressSize, _numBytesInBuffer);

	// newScanBuffer should now contain uncompressed data.
	//delete old buffer, point it at new buffer.
	_saveDataStr.append((char *)newScanBuffer, numDecompressChars);

	delete [] newScanBuffer;
}

void InputStreamMgr::reset()
{
	if (_isBam) {
		return;
	}

	if (!_isStdin) {
		_oldInputStream = _finalInputStream;
		_finalInputStream = new ifstream(_filename.c_str());
	} else {
		if (_isGzipped) {
			_pushBackStreamBuf->pushBack(_compressedSaveData);
		} else {
			_pushBackStreamBuf->pushBack(BTlist<int>(_saveDataStr));
		}
//		_finalInputStream = new istream(_pushBackStreamBuf);
	}
	if (_isGzipped) {
		//the file is gzipped, but is not BAM.
		_infStreamBuf = new InflateStreamBuf(_finalInputStream);
		_oldInputStream = _finalInputStream;
		_finalInputStream = new istream(_infStreamBuf);
	}
}
