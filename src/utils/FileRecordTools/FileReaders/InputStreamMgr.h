/*
 * InputStreamMgr.h
 *
 *  Created on: Mar 21, 2013
 *      Author: nek3d
 */

#ifndef INPUTSTREAMMGR_H_
#define INPUTSTREAMMGR_H_

using namespace std;

#include "PushBackStreamBuf.h"
#include "InflateStreamBuf.h"
#include "QuickString.h"

#include <iostream>

class InputStreamMgr {
public:
	 //Contsructor: 1st arg can be "-" for stdin. set 2nd arg false if fileType already known.
	InputStreamMgr(const QuickString &filename, bool buildScanBuffer = true);
	~InputStreamMgr();
	bool init();

	//use getScanBuffer for auto-detection of file types.
	istream *getFinalStream() { return _finalInputStream; }
	const BTlist<int> &getScanBuffer() const { return _scanBuffer; }
	int getBufferLength() const { return _numBytesInBuffer; }
	void populateScanBuffer();
	void reset();
	const QuickString &getSavedData() const { return _saveDataStr; }
	bool isGzipped() const { return _isGzipped; }
	PushBackStreamBuf *getPushBackStreamBuf() const {return _pushBackStreamBuf; }
	void getSavedData(QuickString &str) const { str = _saveDataStr; }
	bool isBam() const { return _isBam; }

private:
	QuickString _filename;
	PushBackStreamBuf *_pushBackStreamBuf;
	ifstream *_inputFileStream;
	BTlist<int> _scanBuffer;
	QuickString _saveDataStr;
	BTlist<int> _compressedSaveData;
	InflateStreamBuf *_infStreamBuf;
	istream * _finalInputStream;
	istream *_oldInputStream;
	bool _isStdin;
	bool _isGzipped;
	bool _isBam;
	vector<int> _possibleBamCode;
	static const int SCAN_BUFFER_SIZE = 4096; // 4 K buffer
	static const int MIN_SCAN_BUFFER_SIZE = 2048;
	int _numBytesInBuffer; //this will hold the length of the buffer after the scan.

	static const char *FIFO_STRING_LITERAL;
	bool bamDetected(int numChars, int currChar);
	void decompressBuffer();

};

#endif /* INPUTSTREAMMGR_H_ */
