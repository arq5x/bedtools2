
#include "FileRecordMgr.h"
#include "Context.h"
#include "FreeList.h"
#include "Record.h"

FileRecordMgr::FileRecordMgr(int fileIdx, Context *context)
:
  _bufStreamMgr(NULL),
  _contextFileIdx(fileIdx),
  _context(context),
  _fileReader(NULL),
  _fileType(FileRecordTypeChecker::UNKNOWN_FILE_TYPE),
  _recordType(FileRecordTypeChecker::UNKNOWN_RECORD_TYPE),
  _recordMgr(NULL),
  _freeListBlockSize(512),
  _useFullBamTags(false),
  _headerSet(false),
  _prevStart(INT_MAX),
  _prevChromId(-1),
  _mustBeForward(false),
  _mustBeReverse(false),
  _totalRecordLength(0),
  _totalMergedRecordLength(0)
{
}

FileRecordMgr::~FileRecordMgr(){

	if (_bufStreamMgr != NULL) {
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
	}
	if (_fileReader != NULL) {
		close(); //just make sure file was closed.
		delete _fileReader;
		_fileReader = NULL;
	}
	if (_recordMgr != NULL) {
		delete _recordMgr;
		_recordMgr = NULL;
	}
}

bool FileRecordMgr::open(){

	_filename = _context->getInputFileName(_contextFileIdx);
	_bufStreamMgr = new BufferedStreamMgr(_filename);
	if (!_bufStreamMgr->init()) {
		cerr << "Error: unable to open file or unable to determine types for file " << _filename << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}

	_fileType = _bufStreamMgr->getTypeChecker().getFileType();
	_recordType = _bufStreamMgr->getTypeChecker().getRecordType();
	if (_fileType == FileRecordTypeChecker::UNKNOWN_FILE_TYPE || _recordType == FileRecordTypeChecker::UNKNOWN_RECORD_TYPE) {
		cerr << "Error: Unable to determine type for file " << _filename << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}
	allocateFileReader();
	_recordMgr = new RecordMgr(_recordType, _freeListBlockSize);

	_fileReader->setFileName(_filename.c_str());
	_fileReader->setInputStream(_bufStreamMgr);
	_fileReader->setContext(_context);
	if (!_fileReader->open()) {
		cerr << "Error: Types determined but can't open file " << _filename << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}
	_context->setInputFileType(_contextFileIdx, _fileType);
	_context->setInputRecordType(_contextFileIdx, _recordType);

	//if this is a database file, and the numFields is greater than the numFields in other DB files, update.
	// TBD: This is one of many things that will need to be changed to support multiple database files in future versions.
	int numFields = _bufStreamMgr->getTypeChecker().getNumFields();
	if (_contextFileIdx == _context->getDatabaseFileIdx() && numFields > _context->getMaxNumDatabaseFields()) {
		_context->setMaxNumDatabaseFields(numFields);
	}
	if (_fileType == FileRecordTypeChecker::BAM_FILE_TYPE) {
		_context->setHeader(_contextFileIdx, _fileReader->getHeader());
		_context->setReferences(_contextFileIdx, static_cast<BamFileReader *>(_fileReader)->getReferences());
		_headerSet = true;
	}
	return true;
}

void FileRecordMgr::close(){
	if (_bufStreamMgr != NULL) {
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
	}
	if (_fileReader != NULL) {
		_fileReader->close();
		delete _fileReader;
		_fileReader = NULL;
	}
}

bool FileRecordMgr::eof(){
	return _fileReader->eof();
//	return _storedRecords.empty() && _fileReader->eof() ? true:  false;
}

Record *FileRecordMgr::allocateAndGetNextRecord()
{
	if (!_fileReader->isOpen()) {
		return NULL;
	}
	if (!_fileReader->readEntry()) {
		return NULL;
	}
	if (!_headerSet && _fileReader->hasHeader()) {
		_context->setHeader(_contextFileIdx, _fileReader->getHeader());
		_headerSet = true;
	}
	Record *record = NULL;
	record = _recordMgr->allocateRecord();
	if (!record->initFromFile(_fileReader)) {
		_recordMgr->deleteRecord(record);
		return NULL;
	}
	// In the rare case of Bam records where both the read and it's mate failed to map,
	// Ignore the record. Delete and return null
	if (!(record->isUnmapped() && record->isMateUnmapped())) {
		if (!record->coordsValid()) {
			cerr << "Error: Invalid record in file " << _filename << ". Record is " << endl << *record << endl;
			exit(1);
		}

		//test for sorted order, if necessary.
		if (_context->getSortedInput()) {
			testInputSortOrder(record);
		}
	}
	assignChromId(record);
	_totalRecordLength += (unsigned long)(record->getEndPos() - record->getStartPos());
	return record;
}

void FileRecordMgr::assignChromId(Record *record) {
	const QuickString &currChrom = record->getChrName();
	if (currChrom != _prevChrom  && _context->hasGenomeFile()) {
		_prevChromId = _context->getGenomeFile()->getChromId(currChrom);
		record->setChromId(_prevChromId);
	} else {
		record->setChromId(_prevChromId);
	}
}

void FileRecordMgr::testInputSortOrder(Record *record)
{
	// User specified that file must be sorted. Check that it is so.
	// TBD: In future versions, we might not want/need all files to be sorted,
	// even if the -sorted option is used, depending on number of input files
	// and program being run. Should that occur, this block will need adjusting.
	// NEK - 9/5/13


	// Special: For BAM records that aren't mapped, we actually don't want
	// to test the sort order. Another ugly hack sponsored by the letters B, A, and M.
	if (record->isUnmapped()) {
		return;
	}


	const QuickString &currChrom = record->getChrName();
	int currStart = record->getStartPos();
	if (record->isZeroLength()) {
		currStart++;
	}
	if (currChrom != _prevChrom) {
		if ( _foundChroms.find(currChrom) != _foundChroms.end()) {
			//this is a different chrom than the last record had, but we've already seen this chrom.
			sortError(record, false);
		} else {
			//new chrom has not been seen before.
			//TBD: test genome file for ChromId.
			if (_context->hasGenomeFile()) {
				int currChromId = _context->getGenomeFile()->getChromId(currChrom);
				if (currChromId < _prevChromId) {
					sortError(record, true);
				} else {
					_prevChromId = currChromId;
				}
			}
			_foundChroms.insert(currChrom);
			_prevChrom = currChrom;
			_prevStart = INT_MAX;
			record->setChromId(_prevChromId);
		}
	} else if (currStart < _prevStart) { //same chrom as last record, but with lower startPos, so still out of order.
		sortError(record, false);
	}
	_prevStart = currStart;

}

void FileRecordMgr::sortError(const Record *record, bool genomeFileError)
{
	if (genomeFileError) {
		cerr << "Error: Sorted input specified, but the file " << _filename << " has the following record with a different sort order than the genomeFile " <<
				_context->getGenomeFile()->getGenomeFileName() << endl;
	} else {
		cerr << "Error: Sorted input specified, but the file " << _filename << " has the following out of order record" << endl;
	}
	cerr << *record << endl;
	exit(1);
}


void FileRecordMgr::deleteRecord(const Record *record) {
	_recordMgr->deleteRecord(record);
}

void FileRecordMgr::allocateFileReader()
{
	switch (_fileType) {
	case FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE:
	case FileRecordTypeChecker::VCF_FILE_TYPE:
		_fileReader = new SingleLineDelimTextFileReader(_bufStreamMgr->getTypeChecker().getNumFields(), _bufStreamMgr->getTypeChecker().getDelimChar());
		break;

	case FileRecordTypeChecker::BAM_FILE_TYPE:
		_fileReader = new BamFileReader();
		(static_cast<BamFileReader *>(_fileReader))->setUseTags(_context->getUseFullBamTags());
		(static_cast<BamFileReader *>(_fileReader))->setBamReader(_bufStreamMgr->getBamReader());
		break;
	default:
		break;
	}
}


#ifdef false

Record *FileRecordMgr::allocateAndGetNextMergedRecord(WANT_STRAND_TYPE desiredStrand, int maxDistance) {
	RecordKeyList recList;
	if (!allocateAndGetNextMergedRecord(recList, desiredStrand, maxDistance)) {
		return NULL;
	}
	deleteAllMergedItemsButKey(recList);
	return const_cast<Record *>(recList.getKey()); //want key to be non-const
}

bool FileRecordMgr::allocateAndGetNextMergedRecord(RecordKeyList & recList, WANT_STRAND_TYPE desiredStrand, int maxDistance)
{
	if (!recList.allClear()) {
		deleteMergedRecord(recList);
	}

	_mustBeForward = desiredStrand == SAME_STRAND_FORWARD;
	_mustBeReverse = desiredStrand == SAME_STRAND_REVERSE;

	Record *startRecord = tryToTakeFromStorage();

	// if we couldn't use a previously stored record for starters,
	//then begin with a new one that matches strand criteria.
	while (startRecord == NULL) {
		startRecord = allocateAndGetNextRecord();
		if (startRecord == NULL) { //hit EOF!!
			return false;
		}

		if (_mustBeForward && !startRecord->getStrand()) {
			//record is reverse, wanted forward.
			addToStorage(startRecord);
			startRecord = NULL;
		} else if (_mustBeReverse && startRecord->getStrand()) {
			//record is forward, wanted reverse
			addToStorage(startRecord);
			startRecord = NULL;
		}
	}

	// OK!! We have a start record!

	_mustBeForward = desiredStrand == SAME_STRAND_FORWARD || (desiredStrand == SAME_STRAND_EITHER && startRecord->getStrand());
	_mustBeReverse = desiredStrand == SAME_STRAND_REVERSE || (desiredStrand == SAME_STRAND_EITHER && !startRecord->getStrand());

	const QuickString &currChrom = startRecord->getChrName();
	_foundChroms.insert(currChrom);

	bool madeComposite = false;
	recList.push_back(startRecord);
	recList.setKey(startRecord); //key of recList will just be the startRecord unless we're able to merge more.

	bool currStrand = startRecord->getStrand();
	bool mustMatchStrand = desiredStrand != ANY_STRAND;

	int currEnd = startRecord->getEndPos();
	//now look for more records to merge with this one.
	//stop when they're out of range, not on the same chromosome, or we hit EOF.
	//ignore if they don't comply with strand.
	Record *nextRecord = NULL;
	while (nextRecord == NULL) {
		bool takenFromStorage = false;
		nextRecord = mustMatchStrand ? tryToTakeFromStorage(currStrand) : tryToTakeFromStorage();
		if (nextRecord == NULL) {
			nextRecord = allocateAndGetNextRecord();
		} else {
			takenFromStorage = true;
		}
		if (nextRecord == NULL) { // EOF hit
			break;
		}
		const QuickString &newChrom = nextRecord->getChrName();
		if (newChrom != currChrom) { //hit a different chromosome.
			if (_foundChroms.find(newChrom) == _foundChroms.end() || takenFromStorage) {
				//haven't seen this chromosome before.
				addToStorage(nextRecord);
				break;
			} else {
				//different strand, but we've already seen this chrom. File is not sorted.
				fprintf(stderr, "ERROR: Input file %s is not sorted by chromosome, startPos.\n", _context->getInputFileName(_contextFileIdx).c_str());
				deleteRecord(nextRecord);
				deleteMergedRecord(recList);
				exit(1);
			}
		}
		int nextStart = nextRecord->getStartPos();
		//is the record out of range?
		if (nextStart > currEnd + maxDistance) {
			//yes, it's out of range.
			addToStorage(nextRecord);
			break;
		}

		//ok, they're on the same chrom and in range. Are we happy with the strand?
		if (mustMatchStrand && nextRecord->getStrand() != currStrand) {
			//no, we're not.
			addToStorage(nextRecord);
			nextRecord = NULL;
			continue;
		}
		//everything's good! do a merge.
		recList.push_back(nextRecord);
		madeComposite = true;
		int nextEnd = nextRecord->getEndPos();
		if (nextEnd > currEnd) {
			currEnd = nextEnd;
		}
		nextRecord = NULL;
	}
	if (madeComposite) {
		Record *newKey = _recordMgr->allocateRecord();
		(*newKey) = (*startRecord);
		newKey->setEndPos(currEnd);
		recList.setKey(newKey);
	}
	_totalMergedRecordLength += (unsigned long)(recList.getKey()->getEndPos() - recList.getKey()->getStartPos());
	return true;
}

void FileRecordMgr::addToStorage(Record *record) {
	_storedRecords.push(record);
}

Record *FileRecordMgr::tryToTakeFromStorage() {
	Record *record = _storedRecords.empty() ? NULL : const_cast<Record *>(_storedRecords.top());
	if (record != NULL) {
		_storedRecords.pop();
	}
	return record;
}

Record *FileRecordMgr::tryToTakeFromStorage(bool strand) {
	Record *record = NULL;
	if(strand) {
		if (_storedRecords.emptyForward()) {
			return NULL;
		} else {
			record = const_cast<Record *>(_storedRecords.topForward());
			_storedRecords.popForward();
			return record;
		}
	} else {
		if (_storedRecords.emptyReverse()) {
			return NULL;
		} else {
			record = const_cast<Record *>(_storedRecords.topReverse());
			_storedRecords.popReverse();
			return record;
		}
	}
}

void FileRecordMgr::deleteMergedRecord(RecordKeyList &recList)
{
	deleteAllMergedItemsButKey(recList);
	deleteRecord(recList.getKey());
	recList.setKey(NULL);
}

void FileRecordMgr::deleteAllMergedItemsButKey(RecordKeyList &recList) {
	//if the key is also in the list, this method won't delete it.
	for (RecordKeyList::const_iterator_type iter = recList.begin(); iter != recList.end(); iter = recList.next()) {
		if (iter->value() == recList.getKey()) {
			continue;
		}
		deleteRecord(iter->value());
	}
	recList.clearList();
}
#endif
