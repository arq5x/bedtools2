/*
 * InputStreamMgr.h
 *
 *  Created on: Mar 21, 2013
 *      Author: nek3d
 */

#ifndef INPUTSTREAMMGR_H_
#define INPUTSTREAMMGR_H_

#include "PushBackStreamBuf.h"
#include "InflateStreamBuf.h"
#include "string.h"
#include "api/BamReader.h"

#include <htslib/bgzf.h>
#include <iostream>

using namespace std;

class InputStreamMgr {
public:
	 //Contsructor: 1st arg can be "-" for stdin. set 2nd arg false if fileType already known.
	InputStreamMgr(const string &filename, bool buildScanBuffer = true);
	~InputStreamMgr();
	bool init();
	int read(char *data, size_t dataSize);

	//use getScanBuffer for auto-detection of file types.
//	istream *getFinalStream() { return _finalInputStream; }
	const BTlist<int> &getScanBuffer() const { return _scanBuffer; }
	int getBufferLength() const { return _numBytesInBuffer; }
	bool populateScanBuffer();
	const string &getSavedData() const { return _saveDataStr; }
	bool isGzipped() const { return _isGzipped; }
	bool isBGzipped() const { return _isBgzipped; }
	bool isBam() const { return _isBam; }
	bool isCram() const { return _isCram; }

	bool isCompressed() const { return _isGzipped || _isBgzipped || _isBam; }
	PushBackStreamBuf *getPushBackStreamBuf() const {return _pushBackStreamBuf; }
//	void getSavedData(string &str) const { str = _saveDataStr; }
	BamTools::BamReader *getBamReader() { return _bamReader; }
	bool resetStream();
	bool getEofHit() { return _eofHit; }

private:
	string _filename;
	PushBackStreamBuf *_pushBackStreamBuf;
	ifstream *_inputFileStream;
	BTlist<int> _scanBuffer;
	string _saveDataStr;
	InflateStreamBuf *_infStreamBuf;
	istream * _finalInputStream;
	istream *_oldInputStream;
	bool _isStdin;
	bool _isGzipped;
	bool _isBam;
	bool _isCram;
	bool _isBgzipped;
	char *_tmpZipBuf;
	bool _bamRuledOut;
	bool _streamFinished;
	vector<int> _possibleBamCode;
	static const int SCAN_BUFFER_SIZE = 4096; // 4 K buffer
	static const int BAM_SCAN_BUFFER_SIZE = 32768; // 32K
	static const int MIN_SCAN_BUFFER_SIZE = 2048;
	int _numBytesInBuffer; //this will hold the length of the buffer after the scan.
	BamTools::BamReader *_bamReader;
	BGZF* _bgStream;
	bool _eofHit;

	static const char *FIFO_STRING_LITERAL;
	bool readZipChunk();
	bool detectBamOrBgzip(int &numChars, int currChar);
//	void decompressBuffer();

};

#endif /* INPUTSTREAMMGR_H_ */
