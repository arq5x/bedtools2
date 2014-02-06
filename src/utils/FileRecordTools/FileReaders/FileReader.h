#ifndef FILEREADER_H_
#define FILEREADER_H_

using namespace std;

#include <string>
#include <fstream>
#include <map>
#include "QuickString.h"
//#include "ContextBase.h"

#include "BufferedStreamMgr.h"

class FileReader {
public:
	FileReader();
	virtual ~FileReader();
	void setFileName(const string &filename) { _filename = filename; }
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
	virtual const QuickString &getHeader() const =0;
	virtual int getNumFields() const = 0;
protected:
	string _filename;
	BufferedStreamMgr *_bufStreamMgr;

	bool _isFileOpen;
	int _currChromId;

};

#endif // FILEREADER_H_
