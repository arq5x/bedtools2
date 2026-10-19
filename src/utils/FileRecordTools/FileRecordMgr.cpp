
#include "FileRecordMgr.h"
#include "FreeList.h"
#include "Record.h"
#include "NewGenomeFile.h"
#include <typeinfo>

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
  _lineReader(NULL)
 {
}

FileRecordMgr::~FileRecordMgr()
{
	close(); 
	delete _lineReader;
	_lineReader = NULL;
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

bool FileRecordMgr::linesLoadable() const
{
	// a subclass may read records its own way
	if (typeid(*this) != typeid(FileRecordMgr) || _fileReader == NULL || !_fileReader->isOpen()) {
		return false;
	}
	// sorted input and genome files give records chromosome ids that depend
	// on the records read before
	if (_fileType != FileRecordTypeChecker::SINGLE_LINE_DELIM_TEXT_FILE_TYPE || _isSortedInput || _hasGenomeFile) {
		return false;
	}
	switch (_recordType) {
	case FileRecordTypeChecker::BED3_RECORD_TYPE:
	case FileRecordTypeChecker::BED4_RECORD_TYPE:
	case FileRecordTypeChecker::BED5_RECORD_TYPE:
	case FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE:
	case FileRecordTypeChecker::BED6_RECORD_TYPE:
	case FileRecordTypeChecker::BED_PLUS_RECORD_TYPE:
	case FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE:
		return true;
	default:
		return false;
	}
}

// Each class linesLoadable allows prints its fields in file order joined by
// tabs, start and end as signed decimal and the rest as read (BedPlus prints
// a tab after its fixed fields, so it needs a field beyond them); a line has
// exactly the file's number of fields, or readEntry exits.
bool FileRecordMgr::linesPrintAsText() const
{
	const SingleLineDelimTextFileReader *reader = static_cast<const SingleLineDelimTextFileReader *>(_fileReader);
	if (reader->getDelimChar() != '\t') {
		return false;
	}
	int numFields = reader->getNumFields();
	switch (_recordType) {
	case FileRecordTypeChecker::BED3_RECORD_TYPE:
		return numFields == 3;
	case FileRecordTypeChecker::BED4_RECORD_TYPE:
	case FileRecordTypeChecker::BEDGRAPH_RECORD_TYPE:
		return numFields == 4;
	case FileRecordTypeChecker::BED5_RECORD_TYPE:
		return numFields == 5;
	case FileRecordTypeChecker::BED6_RECORD_TYPE:
		return numFields == 6;
	case FileRecordTypeChecker::BED_PLUS_RECORD_TYPE:
		return numFields > 3;
	case FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE:
		return numFields > 6;
	default:
		return false;
	}
}

// Field i of the reader's line as str2chrPos reads it: in place when it is 1
// to 18 digits, else copied into field and parsed, which exits on an illegal
// number. Returns whether the value prints back as the field: 1 to 18 digits
// without a leading zero.
static inline bool lineCoord(const SingleLineDelimTextFileReader *reader, int i, string &field, CHRPOS &val)
{
	const char *s = reader->fieldStart(i);
	size_t n = reader->fieldLength(i);
	if (str2chrPosDigits(s, n, val)) {
		return s[0] != '0' || n == 1;
	}
	reader->getField(i, field);
	val = str2chrPos(field.c_str(), field.size());
	return false;
}

// getNextRecord for a file linesLoadable allows, without the Record. Of the
// classes it allows, only Bed3Interval::initFromFile can fail (str2chrPos on
// columns 2 and 3); the error after it prints the Record, so that one is
// built here too.
bool FileRecordMgr::getNextLine(string &chrName, CHRPOS &start, CHRPOS &end, unsigned *flags)
{
	if (!_fileReader->isOpen()) {
		return false;
	}
	if (!_fileReader->readEntry()) {
		return false;
	}
	SingleLineDelimTextFileReader *reader = static_cast<SingleLineDelimTextFileReader *>(_fileReader);
	if (_lineReader == NULL) {
		_lineReader = new SingleLineDelimTextFileReader(reader->getNumFields(), reader->getDelimChar());
		_lineReader->setFileName(_filename);
		_lineReader->setFileIdx(reader->getFileIdx());
	}
	reader->getField(0, chrName);
	bool canonical = lineCoord(reader, 1, _lineField, start);
	canonical = lineCoord(reader, 2, _lineField, end) && canonical;
	// Record::coordsValid and adjustZeroLength
	if (start < 0 || end < 0 || end < start) {
		Record *record = _recordMgr->allocateRecord();
		record->initFromFile(_fileReader);
		cerr << "Error: Invalid record in file " << _filename << ". Record is " << endl << *record << endl;
		exit(1);
	}
	if (flags != NULL) {
		*flags = (start == end ? LINE_ZERO_LENGTH : 0)
		         | (canonical && !reader->lineHasNul() ? LINE_CANONICAL : 0);
	}
	if (start == end) {
		start--;
		end++;
	}
	// assignChromId leaves the id at _prevChromId, which only sorted input or
	// a genome file changes; recordFromLine sets it
	_totalRecordLength += (unsigned long)(end - start);
	return true;
}

bool FileRecordMgr::bulkLoadable() const
{
	if (!linesLoadable()) {
		return false;
	}
	const SingleLineDelimTextFileReader *reader = static_cast<const SingleLineDelimTextFileReader *>(_fileReader);
	char delim = reader->getDelimChar();
	return reader->linesRead() == 0 && reader->getNumFields() >= 3 && !_bufStreamMgr->isCompressed()
	       && delim != '\n' && delim != '\r' && delim != '\0';
}

void FileRecordMgr::makeLineReader()
{
	if (_lineReader == NULL) {
		SingleLineDelimTextFileReader *reader = static_cast<SingleLineDelimTextFileReader *>(_fileReader);
		_lineReader = new SingleLineDelimTextFileReader(reader->getNumFields(), reader->getDelimChar());
		_lineReader->setFileName(_filename);
		_lineReader->setFileIdx(reader->getFileIdx());
	}
}

// getNextLine for a line given: what follows readEntry there.
bool FileRecordMgr::takeNextLine(const char *line, size_t len, string &chrName, CHRPOS &start, CHRPOS &end,
                                 unsigned *flags)
{
	SingleLineDelimTextFileReader *reader = static_cast<SingleLineDelimTextFileReader *>(_fileReader);
	if (!reader->takeLine(line, len)) {
		return false;
	}
	makeLineReader();
	reader->getField(0, chrName);
	bool canonical = lineCoord(reader, 1, _lineField, start);
	canonical = lineCoord(reader, 2, _lineField, end) && canonical;
	if (start < 0 || end < 0 || end < start) {
		Record *record = _recordMgr->allocateRecord();
		record->initFromFile(_fileReader);
		cerr << "Error: Invalid record in file " << _filename << ". Record is " << endl << *record << endl;
		exit(1);
	}
	if (flags != NULL) {
		*flags = (start == end ? LINE_ZERO_LENGTH : 0)
		         | (canonical && !reader->lineHasNul() ? LINE_CANONICAL : 0);
	}
	if (start == end) {
		start--;
		end++;
	}
	_totalRecordLength += (unsigned long)(end - start);
	return true;
}

void FileRecordMgr::skipLines(size_t n, unsigned long totalLength)
{
	makeLineReader();
	static_cast<SingleLineDelimTextFileReader *>(_fileReader)->skipLines(n);
	_totalRecordLength += totalLength;
}

Record *FileRecordMgr::recordFromLine(const char *line, size_t len)
{
	_lineReader->setLine(line, len);
	Record *record = _recordMgr->allocateRecord();
	record->initFromFile(_lineReader);
	record->coordsValid(); // true: getNextLine checked it; adjusts a zero-length record
	assignChromId(record);
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
	cerr << *record << endl;
	exit(1);
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
		break;

	case FileRecordTypeChecker::BAM_FILE_TYPE:
		_fileReader = new BamFileReader();
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
