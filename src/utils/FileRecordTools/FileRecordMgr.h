/*
 * FileRecordMgr.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef FILERECORDMGR_H_
#define FILERECORDMGR_H_

#include <string>
#include "string.h"
#include <set>
//#include "DualQueue.h"

//include headers for all FileReader and derivative classes.
#include "BufferedStreamMgr.h"
#include "FileReader.h"
#include "SingleLineDelimTextFileReader.h"
#include "BamFileReader.h"

//record manager and all record classes
#include "RecordMgr.h"

#include "RecordKeyVector.h"
#include "BlockMgr.h"

class Record;
class NewGenomeFile;

class FileRecordMgr {
public:
	FileRecordMgr(const string & filename);
	virtual ~FileRecordMgr();
	bool open(bool inheader=false);
	void close();
	virtual bool eof();
	void setFileIdx(int fileIdx) { _fileIdx = fileIdx; }
	int getFileIdx() const { return _fileIdx; }


	//This is an all-in-one method to give the user a new record that is initialized with
	//the next entry in the data file.
	//NOTE!! User MUST pass back the returned pointer to deleteRecord method for cleanup!
	//Also Note! User must check for NULL returned, meaning we failed to get the next record.
	virtual Record *getNextRecord(RecordKeyVector *keyList = NULL);
	void deleteRecord(const Record *);
	virtual void deleteRecord(RecordKeyVector *keyList);

	// Loading without Records, for BinTree. linesLoadable says whether this
	// file allows it: a plain FileRecordMgr on a BED3/4/5/6, bedGraph, BED+
	// or BED6+ text file, read unsorted and without a genome file, so every
	// record gets the same chromosome id. getNextLine reads the next entry
	// as getNextRecord does, with the same errors at the same line, but
	// returns only the chromosome and the start and end after the zero-length
	// adjustment; the line stays in currentLine() until the next read.
	// recordFromLine then builds the Record getNextRecord would have
	// returned for a copy of that line. Given flags, getNextLine also sets
	// LINE_ZERO_LENGTH when the line has start == end, and LINE_CANONICAL
	// when its start and end are canonical decimal (digits, no sign or
	// leading zero, at most 18, so they print back as written) and it holds
	// no NUL. When linesPrintAsText holds too, the Record of a canonical
	// line, with the zero-length adjustment undone, prints as the line.
	enum { LINE_ZERO_LENGTH = 1, LINE_CANONICAL = 2 };
	bool linesLoadable() const;
	bool linesPrintAsText() const;
	bool getNextLine(string &chrName, CHRPOS &start, CHRPOS &end, unsigned *flags = NULL);
	// the file index of the Records recordFromLine builds
	int lineFileIdx() const { return _fileReader->getFileIdx(); }
	const char *currentLine() const {
		return static_cast<const SingleLineDelimTextFileReader *>(_fileReader)->currentLine();
	}
	size_t currentLineLength() const {
		return static_cast<const SingleLineDelimTextFileReader *>(_fileReader)->currentLineLength();
	}
	Record *recordFromLine(const char *line, size_t len);
	// The reader of a text file; its current line is the entry last read.
	const SingleLineDelimTextFileReader *textReader() const {
		return static_cast<const SingleLineDelimTextFileReader *>(_fileReader);
	}

	// Loading lines split from the file elsewhere (BinTree's parallel load),
	// for a file linesLoadable allows. bulkLoadable says whether this one
	// can: nothing read from it yet, a delimiter getLine's line scan takes,
	// and no decompression, whose errors print while the file is read, so
	// that reading ahead would print them early. readRaw then gives the
	// file's bytes from its start
	// (BufferedStreamMgr::readRaw), and each line split from them goes,
	// in order, either to takeNextLine, which reads it as getNextLine would
	// have read it from the file (SingleLineDelimTextFileReader::takeLine),
	// or to skipLines, for lines found to be entries elsewhere: n of them,
	// whose end - start, after the zero-length adjustment, sum to
	// totalLength.
	bool bulkLoadable() const;
	size_t readRaw(char *dst, size_t cap) { return _bufStreamMgr->readRaw(dst, cap); }
	bool takeNextLine(const char *line, size_t len, string &chrName, CHRPOS &start, CHRPOS &end,
	                  unsigned *flags = NULL);
	void skipLines(size_t n, unsigned long totalLength);



	const string &getFileName() const { return _filename;}
	bool hasHeader() const { return _fileReader->hasHeader(); }
	const string &getHeader() const { return _fileReader->getHeader(); }

	bool recordsHaveName() const {
		return _bufStreamMgr->getTypeChecker().recordTypeHasName(_recordType);
	}
	bool recordsHaveScore() const {
		return _bufStreamMgr->getTypeChecker().recordTypeHasScore(_recordType);
	}
	bool recordsHaveStrand() const {
		return _bufStreamMgr->getTypeChecker().recordTypeHasStrand(_recordType);
	}

	FileRecordTypeChecker::FILE_TYPE getFileType() const {
		return _fileType;
	}
	FileRecordTypeChecker::RECORD_TYPE getRecordType() const {
		return _recordType;
	}
	const string &getFileTypeName() const {
		return _bufStreamMgr->getTypeChecker().getFileTypeName();
	}

	const string &getRecordTypeName() const {
		return _bufStreamMgr->getTypeChecker().getRecordTypeName();
	}

	const BamTools::RefVector &getBamReferences();
	refs_t* getCramRefs();

	int getNumFields() const { return _fileReader->getNumFields(); }

	//File statistics
	unsigned long getTotalRecordLength() const { return _totalRecordLength; } //sum of length of all returned records
	unsigned long getTotalMergedRecordLength() const { return _totalMergedRecordLength; } // sum of all merged intervals



	//Setting the freeListBlockSize is optional. If the user never calls this,
	//the blockSize defaults to 512.
	void setFreeListBlockSize(int blockSize) { _freeListBlockSize = blockSize; }

	//special: For BAM files, our default is to not use all the
	//tag information in a BAM file, which reduces the run time in some
	//cases by more than 50%. But setting this method to true
	//will use all the tags, if more information is desired.
	//MUST BE CALLED BEFORE THE BAM FILE IS OPEN.
	void setFullBamFlags(bool flag) { _useFullBamTags = flag; }
	void setGenomeFile(NewGenomeFile *genomeFile) {
		_genomeFile = genomeFile;
		_hasGenomeFile = true;
	}

	void setIsSorted(bool val) { _isSortedInput = val; }
	void setIoBufSize(int val) { _ioBufSize = val; }
	void setNoEnforceCoordSort(bool val) { _noEnforceCoordSort = val; }
	void setIsGroupBy(bool val) { _isGroupBy = val; }

	bool isCram() const { return _isCram; }

protected:
	int _fileIdx;
	string _filename;
	BufferedStreamMgr *_bufStreamMgr;

	FileReader *_fileReader;
	FileRecordTypeChecker::FILE_TYPE _fileType;
	FileRecordTypeChecker::RECORD_TYPE _recordType;
	RecordMgr *_recordMgr;
	bool _isSortedInput;
	int _freeListBlockSize;
	bool _useFullBamTags;

	//members for enforcing sorted order.
        std::set<string> _foundChroms;
	string _prevChrom;
	CHRPOS _prevStart;
	int _prevChromId;


	bool _mustBeForward;
	bool _mustBeReverse;

	//available stats after run
	unsigned long _totalRecordLength;
	unsigned long _totalMergedRecordLength;

	BlockMgr *_blockMgr;
	BamTools::BamReader *_bamReader;
	bool _hasGenomeFile;
	NewGenomeFile *_genomeFile;
	int _ioBufSize;
	bool _noEnforceCoordSort; //only true for GroupBy
	bool _isGroupBy; //hopefully also only true for GroupBy
	bool _isCram;

	// parses the lines recordFromLine is given; its own copy, so the file's
	// reader is never disturbed
	SingleLineDelimTextFileReader *_lineReader;
	string _lineField;
	void makeLineReader();

	void allocateFileReader(bool inheader=false);
	void testInputSortOrder(Record *record);
	void assignChromId(Record *);
	void sortError(const Record *record, bool genomeFileError);
};


#endif /* FILERECORDMGR_H_ */
