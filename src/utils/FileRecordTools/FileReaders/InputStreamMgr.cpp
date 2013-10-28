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

const char *InputStreamMgr::FIFO_STRING_LITERAL = "/dev/fd";

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
 _isBgzipped(false),
 _bamRuledOut(false),
 _numBytesInBuffer(0),
 _bamReader(NULL),
 _bgStream(NULL)
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
	if (_bamReader != NULL) {
		delete _bamReader;
		_bgStream = NULL;
	}
	if (_bgStream != NULL) {
		delete _bgStream;
		_bgStream = NULL;
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
		if (strncmp(_filename.c_str(), FIFO_STRING_LITERAL, strlen(FIFO_STRING_LITERAL)) == 0) {
			_isStdin = true;
		}
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

int InputStreamMgr::read(char *data, size_t dataSize)
{
	if (_isBgzipped) {
		return (int)(_bgStream->Read(data, dataSize));
	}
	_finalInputStream->read(data, dataSize);
	return _finalInputStream->gcount();
}

void InputStreamMgr::populateScanBuffer()
{
	_scanBuffer.clear();
	int numChars=0;
	int currChar = 0;
	bool mustAppend = true;
	while (1) {
		mustAppend = true;
		currChar = _pushBackStreamBuf->sbumpc();
		//Stop when EOF hit.
		if (currChar == EOF) {
			break;
		}
		numChars++;
		_scanBuffer.push_back(currChar);
		if (_isGzipped) {
			if (!_bamRuledOut && detectBamOrBgzip(numChars, currChar, mustAppend)) {
				return;
			}
			if (numChars == 0) {
				continue; //this will only happen when we've just discovered that this
				//is definitely not BAM, and want to start over.
			}
			if (mustAppend) {
				_compressedSaveData.push_back(currChar);
			}
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

bool InputStreamMgr::detectBamOrBgzip(int &numChars, int currChar, bool &mustAppend)
{
	//Look for the BAM magic string "BAM\1" in the first fouur characters of the input stream.
	//In compressed form, the first char is the gzip signifier, which was already found.
	//The next three are the integers 139, 8, and 4.
	if (numChars < 5) {
		_possibleBamCode[numChars -1] = currChar;
		//special: test for BAM
		if (numChars == 4 && _possibleBamCode[1] == 139 && _possibleBamCode[2] == 8 && _possibleBamCode[3] == 4) {
			//BAM magic string detected.This is either a BAM or bgzip file. To find out which, we have to try and
			//open the file as BAM, with a BAM reader, and see if the header and references are both non-empty.
			//However, if they are empty, we will have had to save all bytes consumed in the attempt, meaning still
			//fill the scanBuffer and push it back onto the pushBackStream as normal.
			for (; numChars < BAM_SCAN_BUFFER_SIZE; numChars++) {
				currChar = _pushBackStreamBuf->sbumpc();
				//Stop when EOF hit.
				if (currChar == EOF) {
					break;
				}
				_scanBuffer.push_back(currChar);

			}
			_pushBackStreamBuf->pushBack(_scanBuffer);

			//ok, now all the data read so far is saved in the scan buffer, and pushbackstream is reset.
			//now we make a BamReader and try to open the file.


			_bamReader = new BamTools::BamReader();
			_bamReader->OpenStream(_finalInputStream);
			QuickString bamHeader(_bamReader->GetHeaderText());
			BamTools::RefVector bamRefs(_bamReader->GetReferenceData());

			if (bamHeader.empty() || bamRefs.empty()) {
				//This is NOT a bam file.
				_pushBackStreamBuf->clear();
				_compressedSaveData.clear();
				//Put all bytes read so far back onto the scan buffer, then reset
				//everything so that we're effectively starting over.
				_pushBackStreamBuf->pushBack(_scanBuffer);
				_scanBuffer.clear();
				mustAppend = false;
				numChars = 0;
				_isBam = false;
				_isBgzipped = true;
				_bamRuledOut = true;
				_numBytesInBuffer = 0;
				delete _bamReader;
				_bamReader = NULL;
				return false;
			}
			_isBam = true;
			_numBytesInBuffer = _scanBuffer.size();
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
		//For file input, just re-open the file.
		_oldInputStream = _finalInputStream;
		_finalInputStream = new ifstream(_filename.c_str());
	} else {
		if (_isBgzipped) {
			for (BTlist<int>::const_iterator_type iter = _pushBackStreamBuf->_buffer.begin();
					iter != _pushBackStreamBuf->_buffer.end(); iter = _pushBackStreamBuf->_buffer.next()) {
				_compressedSaveData.push_back(iter->value());
			}
			_pushBackStreamBuf->clear();
			_pushBackStreamBuf->pushBack(_compressedSaveData);
		} else if (_isGzipped) {
			_pushBackStreamBuf->pushBack(_compressedSaveData);
		} else {
			_pushBackStreamBuf->pushBack(BTlist<int>(_saveDataStr));
		}
//		_finalInputStream = new istream(_pushBackStreamBuf);
	}
	if (_isBgzipped) {
		//The file is bgzipped, but not BAM.
		_bgStream = new BamTools::Internal::BgzfStream();
		_bgStream->OpenStream(_finalInputStream, BamTools::IBamIODevice::ReadOnly);
	} else if (_isGzipped) {
		//the file is gzipped, but is not bgzipped or BAM.
		_infStreamBuf = new InflateStreamBuf(_finalInputStream);
		if (_oldInputStream != NULL) {
			delete _oldInputStream;
		}
		_oldInputStream = _finalInputStream;
		_finalInputStream = new istream(_infStreamBuf);
	}
}


