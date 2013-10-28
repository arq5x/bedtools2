#ifndef FILEREADER_H_
#define FILEREADER_H_

using namespace std;

#include <string>
#include <fstream>
#include <map>
#include "QuickString.h"
#include "Context.h"

#include "BufferedStreamMgr.h"

class FileReader {
public:
	FileReader();
	virtual ~FileReader();
	void setFileName(const string &filename) { _filename = filename; }
//	void setInputStream(istream *inputStream) {
//		_inputStream = inputStream;
//		_externalInputStream = true;
//	}
	void setInputStream(BufferedStreamMgr *bufStreamMgr) {
		_bufStreamMgr = bufStreamMgr;
//		_inputStream = _bufStreamMgr->getStream();
//		_externalInputStream = true;
//		_useBufStream = true;

		// This will short circuit the open method. BufferedStreamMgr does it's own file opening.
		//However, for BAM, we want to re-open it.
		_isFileOpen = true; //_bufStreamMgr->getTypeChecker().isBam() ? false : true;

	}
	void setContext(const Context *context) { _context = context; }
	virtual bool open(); //open the file
	virtual bool isOpen() const { return _isFileOpen; }
	virtual void close();
	virtual bool eof() const;
	virtual bool readEntry() =0; // this is an abstract base class.
	//Derived classes could read text, BAM, possibly FASTA, and/or Variant files, as well as Blast/SIM4 output.
	virtual int getCurrChromdId() const { return _currChromId; }
	virtual bool hasHeader() const = 0;
	virtual const QuickString &getHeader() const =0;
protected:
	string _filename;
//	istream   *_inputStream;
	BufferedStreamMgr *_bufStreamMgr;

	bool _isFileOpen;
//	bool _mustDeleteInputStream;
//	bool _externalInputStream;
//	bool _useBufStream;
	int _currChromId;
	Context const *_context;

};

#endif // FILEREADER_H_
