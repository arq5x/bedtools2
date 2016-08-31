#ifndef FILEREADER_H_
#define FILEREADER_H_

#include <string>
#include <fstream>
#include <map>
#include "string.h"
//#include "ContextBase.h"

#include "BufferedStreamMgr.h"

using namespace std;

class FileReader {
public:
	FileReader();
	virtual ~FileReader();
	void setFileName(const string &filename) { _filename = filename; }
	virtual int getFileIdx() const { return _fileIdx; }
	virtual void setFileIdx(int fileIdx) { _fileIdx = fileIdx; }

	void setInputStream(BufferedStreamMgr *bufStreamMgr) {
		_bufStreamMgr = bufStreamMgr;
		_isFileOpen = true;

	}
	virtual bool open(); //open the file
	virtual bool isOpen() const { return _isFileOpen; }
	virtual void close();
	virtual bool eof() const;
	virtual bool readEntry() =0; // this is an abstract base class.
	virtual int getCurrChromdId() const { return _currChromId; }
	virtual bool hasHeader() const = 0;
	virtual const string &getHeader() const =0;
	virtual int getNumFields() const = 0;
protected:
	int _fileIdx;
	string _filename;
	BufferedStreamMgr *_bufStreamMgr;

	bool _isFileOpen;
	int _currChromId;

};

#endif // FILEREADER_H_
