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

	void allocateFileReader(bool inheader=false);
	void testInputSortOrder(Record *record);
	void assignChromId(Record *);
	void sortError(const Record *record, bool genomeFileError);
};


#endif /* FILERECORDMGR_H_ */
