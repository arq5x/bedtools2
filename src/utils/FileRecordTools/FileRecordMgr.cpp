
#include "FileRecordMgr.h"
#include "FreeList.h"
#include "Record.h"
#include "NewGenomeFile.h"

FileRecordMgr::FileRecordMgr(const string &filename)
: _fileIdx(-1),
  _filename(filename),
  _bufStreamMgr(NULL),
  _fileReader(NULL),
  _fileType(FileRecordTypeChecker::UNKNOWN_FILE_TYPE),
  _recordType(FileRecordTypeChecker::UNKNOWN_RECORD_TYPE),
  _recordMgr(NULL),
  _isSortedInput(false),
  _freeListBlockSize(512),
  _useFullBamTags(false),
  _prevStart(INT_MAX),
  _prevChromId(-1),
  _mustBeForward(false),
  _mustBeReverse(false),
  _totalRecordLength(0),
  _totalMergedRecordLength(0),
  _blockMgr(NULL),
  _bamReader(NULL),
  _hasGenomeFile(false),
  _genomeFile(NULL),
  _ioBufSize(0),
  _noEnforceCoordSort(false),
  _isGroupBy(false),
  _isCram(false),
  _coordsOnly(false),
  _sweepCoords(false),
  _queryFast(false),
  _queryChrSet(false)
 {
}

FileRecordMgr::~FileRecordMgr()
{
	close(); 
	//delete _recordMgr;
	_recordMgr = NULL;
}

bool FileRecordMgr::open(bool inheader){
	_bufStreamMgr = new BufferedStreamMgr(_filename);
	_bufStreamMgr->getTypeChecker().setInHeader(inheader);

	if (_ioBufSize > 0) _bufStreamMgr->setIoBufSize(_ioBufSize);
	
	if (_isGroupBy) {
		_bufStreamMgr->getTypeChecker().setIsGroupBy(true);
	}
	
	if (!_bufStreamMgr->init()) {
		cerr << "Error: unable to open file or unable to determine types for file " << _filename << endl;
		cerr << endl;
		cerr << "- Please ensure that your file is TAB delimited (e.g., cat -t FILE)." << endl;
		cerr << "- Also ensure that your file has integer chromosome coordinates in the " << endl
		     << "  expected columns (e.g., cols 2 and 3 for BED)." << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}
	_fileType = _bufStreamMgr->getTypeChecker().getFileType();
	_recordType = _bufStreamMgr->getTypeChecker().getRecordType();

	//HACK: If groupBy and not Bam, over-ride file type.
	if (_isGroupBy && _fileType != FileRecordTypeChecker::BAM_FILE_TYPE) 
	{
		_bufStreamMgr->getTypeChecker().setFileType(FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE);
		_bufStreamMgr->getTypeChecker().setRecordType(FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE);
		_fileType = FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE;
		_recordType = FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE;
	}

	if(_fileType == FileRecordTypeChecker::BAM_FILE_TYPE)
		_isCram = _bufStreamMgr->getTypeChecker().isCram();

	if (_fileType == FileRecordTypeChecker::UNKNOWN_FILE_TYPE || _recordType == FileRecordTypeChecker::UNKNOWN_RECORD_TYPE) {
		cerr << "Error: Unable to determine type for file " << _filename << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}
	allocateFileReader(inheader);
	_recordMgr = new RecordMgr(_recordType, _freeListBlockSize);

	_fileReader->setFileName(_filename.c_str());
	_fileReader->setInputStream(_bufStreamMgr);
	if (!_fileReader->open()) {
		cerr << "Error: Types determined but can't open file " << _filename << endl;
		delete _bufStreamMgr;
		_bufStreamMgr = NULL;
		exit(1);
	}

	return true;
}

void FileRecordMgr::close(){
	delete _bufStreamMgr;
	_bufStreamMgr = NULL;

	if (_fileReader != NULL) {
		_fileReader->close();
		delete _fileReader;
		_fileReader = NULL;
	}
}

bool FileRecordMgr::eof(){
	return _fileReader->eof();
}

Record *FileRecordMgr::getNextRecord(RecordKeyVector *keyList)
{
	if (!_fileReader->isOpen()) {
		return NULL;
	}
	if (!_fileReader->readEntry()) {
		return NULL;
	}
	return recordFromLine(keyList);
}

Record *FileRecordMgr::getNextQueryRecord()
{
	bool fast = _queryFast && !_hasGenomeFile && !_isSortedInput && !_coordsOnly
		&& (_recordType == FileRecordTypeChecker::BED3_RECORD_TYPE
			|| _recordType == FileRecordTypeChecker::BED6_RECORD_TYPE
			|| _recordType == FileRecordTypeChecker::BED_PLUS_RECORD_TYPE
			|| _recordType == FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE);
	if (fast) {
		int nf = _fileReader->getNumFields();
		fast = _recordType == FileRecordTypeChecker::BED3_RECORD_TYPE ? nf == 3
			: _recordType == FileRecordTypeChecker::BED6_RECORD_TYPE ? nf == 6
			: _recordType == FileRecordTypeChecker::BED_PLUS_RECORD_TYPE ? nf > 3
			: nf > 6;
	}
	if (!fast) {
		return getNextRecord();
	}
	if (!_fileReader->isOpen()) {
		return NULL;
	}
	SingleLineDelimTextFileReader *rr = static_cast<SingleLineDelimTextFileReader*>(_fileReader);
	if (_queryChrSet) {
		CHRPOS s, e;
		int got = rr->readCoordsLine(_queryChr.str(), s, e);
		if (got == 0) {
			return NULL;
		}
		if (got == 1) {
			Record *record = _recordMgr->allocateRecord();
			if (record->initFastLine(_queryChr, s, e, rr->getFileIdx(), rr->getLineData(), rr->getLineLength(),
					rr->getCoordsTail(), rr->getNumFields())) {
				record->setChromId(_prevChromId);
				_totalRecordLength += (unsigned long)(record->getEndPos() - record->getStartPos());
				record->setFileRecordManager(this);
				return record;
			}
			// the gate above makes this unreachable
			cerr << "Error: unable to read a record from file " << _filename << endl;
			exit(1);
		}
	} else if (!_fileReader->readEntry()) {
		return NULL;
	}
	Record *record = recordFromLine(NULL);
	if (record != NULL) {
		_queryChr = record->getChrRef();
		_queryChrSet = true;
	}
	return record;
}

Record *FileRecordMgr::getNextRecordSweep()
{
	if (!_sweepCoords) return getNextRecord();
	if (!_fileReader->isOpen() || !_fileReader->readEntry()) {
		return NULL;
	}
	SingleLineDelimTextFileReader *r = static_cast<SingleLineDelimTextFileReader*>(_fileReader);
	const int *dp = r->getDelimPositions();
	size_t chrLen = (size_t)(dp[1] - dp[0] - 1);
	bool same = chrLen == _sweepChr.size() && memcmp(r->getLineData() + dp[0] + 1, _sweepChr.data(), chrLen) == 0;
	r->setCoordsOnly(same);
	Record *rec = recordFromLine(NULL);
	if (!same && rec != NULL) _sweepChr = rec->getChrName();
	return rec;
}

void FileRecordMgr::getCurrentLine(const char *&line, size_t &lineLen) const
{
	const SingleLineDelimTextFileReader *r = static_cast<SingleLineDelimTextFileReader*>(_fileReader);
	line = r->getLineData();
	lineLen = r->getLineLength();
}

int FileRecordMgr::getNextCoords(const string *lastChr, CHRPOS &start, CHRPOS &end, bool &zeroLen, Record *&rec,
		const char **lineOut, size_t *lineLenOut)
{
	rec = NULL;
	if (!_fileReader->isOpen()) {
		return 0;
	}
	if (lastChr != NULL && !_hasGenomeFile && !_isSortedInput) {
		SingleLineDelimTextFileReader *rr = static_cast<SingleLineDelimTextFileReader*>(_fileReader);
		CHRPOS fs, fe;
		int got = rr->readCoordsLine(*lastChr, fs, fe);
		if (got == 0) {
			return 0;
		}
		if (got == 1) {
			bool zl = (fs == fe);
			if (zl) {
				fs--;
				fe++;
			}
			zeroLen = zl;
			start = fs;
			end = fe;
			_totalRecordLength += (unsigned long)(fe - fs);
			if (lineOut != NULL) {
				*lineOut = rr->getLineData();
				*lineLenOut = rr->getLineLength();
			}
			return 1;
		}
	}
	else if (!_fileReader->readEntry()) {
		return 0;
	}
	if (lastChr != NULL && !_hasGenomeFile && !_isSortedInput) {
		const SingleLineDelimTextFileReader *r = static_cast<SingleLineDelimTextFileReader*>(_fileReader);
		const char *line = r->getLineData();
		const int *dp = r->getDelimPositions();
		size_t chrLen = (size_t)(dp[1] - dp[0] - 1);
		CHRPOS s, e;
		if (chrLen == lastChr->size() && memcmp(line + dp[0] + 1, lastChr->data(), chrLen) == 0
				&& r->getPlainPos(1, s) && r->getPlainPos(2, e) && s <= e) {
			zeroLen = (s == e);
			if (zeroLen) {
				s--;
				e++;
			}
			start = s;
			end = e;
			_totalRecordLength += (unsigned long)(e - s);
			if (lineOut != NULL) {
				*lineOut = line;
				*lineLenOut = r->getLineLength();
			}
			return 1;
		}
	}
	rec = recordFromLine(NULL);
	return rec == NULL ? 0 : 2;
}

Record *FileRecordMgr::recordFromLine(RecordKeyVector *keyList)
{
	Record *record = NULL;
	record = _recordMgr->allocateRecord();
	if (!record->initFromFile(_fileReader)) {
		_recordMgr->deleteRecord(record);
		return NULL;
	}

	// If the record is unmapped, don't test for valid coords or sort order,
	// but still return it so the -v (noHit) option and the like will still
	// see it.

	if (!record->isUnmapped() ) {
		if (!record->coordsValid() && (record->getType() != FileRecordTypeChecker::NO_POS_PLUS_RECORD_TYPE)) {
			fullReparse(record);
			cerr << "Error: Invalid record in file " << _filename << ". Record is " << endl << *record << endl;
			exit(1);
		}

		//test for sorted order, if necessary.
		if (_isSortedInput && !_noEnforceCoordSort) {
			testInputSortOrder(record);
		}
	}
	assignChromId(record);
	_totalRecordLength += (unsigned long)(record->getEndPos() - record->getStartPos());
	if (keyList != NULL) {
		keyList->setKey(record);
	}
	record->setFileRecordManager(this);
	return record;
}

void FileRecordMgr::assignChromId(Record *record) {
	const string &currChrom = record->getChrName();
	if (currChrom != _prevChrom  && _hasGenomeFile) {
		_prevChromId = _genomeFile->getChromId(currChrom);
		record->setChromId(_prevChromId);
	} else {
		record->setChromId(_prevChromId);
	}
}

void FileRecordMgr::testInputSortOrder(Record *record)
{

	// Special: For BAM records that aren't mapped, we actually don't want
	// to test the sort order. Another ugly hack sponsored by the letters B, A, and M.
	if (record->isUnmapped()) {
		return;
	}


	const string &currChrom = record->getChrName();
	CHRPOS currStart = record->getStartPos();
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
			if (_hasGenomeFile) {
				int currChromId = _genomeFile->getChromId(currChrom);
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
				_genomeFile->getGenomeFileName() << endl;
	} else {
		cerr << "Error: Sorted input specified, but the file " << _filename << " has the following out of order record" << endl;
	}
	fullReparse(const_cast<Record *>(record));
	cerr << *record << endl;
	exit(1);
}

// An error message prints the whole record: parse the current line in full first.
void FileRecordMgr::fullReparse(Record *record)
{
	SingleLineDelimTextFileReader *r = dynamic_cast<SingleLineDelimTextFileReader *>(_fileReader);
	if (r == NULL || !r->getCoordsOnly()) return;
	r->setCoordsOnly(false);
	record->initFromFile(r);
}


void FileRecordMgr::deleteRecord(const Record *record) {
	_recordMgr->deleteRecord(record);
}

void FileRecordMgr::deleteRecord(RecordKeyVector *keyList) {
	_recordMgr->deleteRecord(keyList->getKey());
}

void FileRecordMgr::allocateFileReader(bool inheader)
{
	switch (_fileType) {
	case FileRecordTypeChecker::EMPTY_FILE_TYPE:
	case FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE:
	case FileRecordTypeChecker::VCF_FILE_TYPE:
		_fileReader = new SingleLineDelimTextFileReader(_bufStreamMgr->getTypeChecker().getNumFields(), _bufStreamMgr->getTypeChecker().getDelimChar());
		static_cast<SingleLineDelimTextFileReader *>(_fileReader)->setInHeader(inheader);
		static_cast<SingleLineDelimTextFileReader *>(_fileReader)->setCoordsOnly(_coordsOnly);
		break;

	case FileRecordTypeChecker::BAM_FILE_TYPE:
		_fileReader = new BamFileReader();
		_sweepCoords = false;
		(static_cast<BamFileReader *>(_fileReader))->setUseTags(_useFullBamTags);
		(static_cast<BamFileReader *>(_fileReader))->setBamReader(_bufStreamMgr->getBamReader());
		break;
	default:
		break;
	}
	_fileReader->setFileIdx(_fileIdx);
}

const BamTools::RefVector & FileRecordMgr::getBamReferences() {
	// exta safety check to insure user checked the file type first.
	if (_fileType != FileRecordTypeChecker::BAM_FILE_TYPE) {
		cerr << "Error: Attempted to get BAM references from file " << _filename << ", which is NOT a BAM file." << endl;
		exit(1);
	}
	return static_cast<BamFileReader *>(_fileReader)->getReferences();
}
refs_t* FileRecordMgr::getCramRefs() {
	if (_fileType != FileRecordTypeChecker::BAM_FILE_TYPE) {
		cerr << "Error: Attempted to get BAM references from file " << _filename << ", which is NOT a BAM file." << endl;
		exit(1);
	}
	return static_cast<BamFileReader *>(_fileReader)->getCramRefs();
}
