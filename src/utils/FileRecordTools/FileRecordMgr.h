/*
 * FileRecordMgr.h
 *
 *  Created on: Nov 8, 2012
 *      Author: nek3d
 */

#ifndef FILERECORDMGR_H_
#define FILERECORDMGR_H_

using namespace std;

#include <string>
#include "QuickString.h"
#include <set>
//#include "DualQueue.h"

//include headers for all FileReader and derivative classes.
#include "Context.h"
#include "BufferedStreamMgr.h"
#include "FileReader.h"
#include "SingleLineDelimTextFileReader.h"
#include "BamFileReader.h"

//record manager and all record classes
#include "RecordMgr.h"

#include "RecordKeyList.h"
#include "BlockMgr.h"

class Record;
class FileRecordMgr {
public:
	FileRecordMgr(int inputFileIdx, Context *context);
	~FileRecordMgr();
	bool open();
	void close();
	bool eof();

	//NOTE: FOR Non-BAM FILES, HEADER INFORMATION IS CURRENTLY ONLY AVAILABLE AFTER 1st CALL
	//TO allocateAndGetNextRecord().
	bool hasHeader() const { return _fileReader->hasHeader(); }
	const QuickString &getHeader() const { return _fileReader->getHeader(); }

	bool recordsHaveName() const {
		testBufferedMgr();
		return _bufStreamMgr->getTypeChecker().recordTypeHasName(_recordType);
	}
	bool recordsHaveScore() const {
		testBufferedMgr();
		return _bufStreamMgr->getTypeChecker().recordTypeHasScore(_recordType);
	}
	bool recordsHaveStrand() const {
		testBufferedMgr();
		return _bufStreamMgr->getTypeChecker().recordTypeHasStrand(_recordType);
	}

	FileRecordTypeChecker::FILE_TYPE getFileType() const {
		testBufferedMgr();
		return _fileType;
	}
	FileRecordTypeChecker::RECORD_TYPE getRecordType() const {
		testBufferedMgr();
		return _recordType;
	}
	const string &getFileTypeName() const {
		testBufferedMgr();
		return _bufStreamMgr->getTypeChecker().getFileTypeName();
	}

	const string &getRecordTypeName() const {
		testBufferedMgr();
		return _bufStreamMgr->getTypeChecker().getRecordTypeName();
	}

	//This is an all-in-one method to give the user a new record that is initialized with
	//the next entry in the data file.
	//NOTE!! User MUST pass back the returned pointer to deleteRecord method for cleanup!
	//Also Note! User must check for NULL returned, meaning we failed to get the next record.
	Record *allocateAndGetNextRecord();
	void deleteRecord(const Record *);

#ifdef false
	//////////////////////////////////////////////////////////////////////////////////
	//
	// 			MERGED RECORDS
	//
	//this will give a single "meta" record containing "flattened" or merged records.
	//
	// 1st ARG: Pass an empty RecordKeyList. When done, will have a pair: 1st is the final merged record,
	//			second is list of constituent Records merged.
	//			** NOTE ** If the RecordKeyList is not empty, this method will empty it for you and delete all contents!
	//
	// 2nd ARG: Choose from WANT_STRAND_TYPE, defined below below
	//
	// 3rd ARG: allows for nearby records, i.e. maxDistance 100 will merge records <= 100 bases apart. Default 0 means only
	//			merge records that actually intersect.
	//
	// Return value: true if any records found. False if eof hit before records matching requested parameters found.

	typedef enum { SAME_STRAND_FORWARD, //must all be forward strand
			SAME_STRAND_REVERSE, //must all be reverse strand
			SAME_STRAND_EITHER, //must be same strand, but can be either forward or reverse
			ANY_STRAND } //do no care about strand (Default value)
	WANT_STRAND_TYPE;

	//
	// WARNING!! Specifying a strand will keep all records on the other strand in memory!!
	// This is done so that requests for records on that other strand can still be met.
	// For now, use this method at any time to purge the kept records from memory, such as
	// when changing chromosomes, for example.
	void purgeKeepList();
	bool allocateAndGetNextMergedRecord(RecordKeyList & recList, WANT_STRAND_TYPE desiredStrand = ANY_STRAND, int maxDistance = 0);
	void deleteMergedRecord(RecordKeyList &recList); // MUST use this method for cleanup!

	//this method will allocate a new record of merged records, but the returned record should only be passed to the deleteRecord method
	//for cleanup, not to the delete mmerged record.
	Record *allocateAndGetNextMergedRecord(WANT_STRAND_TYPE desiredStrand = ANY_STRAND, int maxDistance = 0);

	//
	// 				END MERGED RECORDS
	//
	//////////////////////////////////////////////////////////////////////////////////
#endif

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

private:
	QuickString _filename;
	BufferedStreamMgr *_bufStreamMgr;
	int _contextFileIdx;

	Context *_context;
	FileReader *_fileReader;
	FileRecordTypeChecker::FILE_TYPE _fileType;
	FileRecordTypeChecker::RECORD_TYPE _recordType;
	RecordMgr *_recordMgr;
	int _freeListBlockSize;
	bool _useFullBamTags;
	bool _headerSet;

	//members for enforcing sorted order.
	set<QuickString> _foundChroms;
	QuickString _prevChrom;
	int _prevStart;
	int _prevChromId;

	//members for handling merged records
//	DualQueue<Record *, DualQueueAscending > _storedRecords;

	bool _mustBeForward;
	bool _mustBeReverse;

	//available stats after run
	unsigned long _totalRecordLength;
	unsigned long _totalMergedRecordLength;

	BlockMgr *_blockMgr;
	BamTools::BamReader *_bamReader;


	void allocateFileReader();
	void testInputSortOrder(Record *record);
	void assignChromId(Record *);
	void sortError(const Record *record, bool genomeFileError);
	void testBufferedMgr() const {
		if (_bufStreamMgr ==  NULL) {
			cerr << "Error: attempted to access type information for file " << _context->getInputFileName(_contextFileIdx) << " before it was opened." << endl;
			exit (1);
		}
	}

#ifdef false
	void deleteAllMergedItemsButKey(RecordKeyList &recList);
	void addToStorage(Record *record);
	Record *tryToTakeFromStorage();
	Record *tryToTakeFromStorage(bool strand);
#endif


};


#endif /* FILERECORDMGR_H_ */
