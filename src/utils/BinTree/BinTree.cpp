#include "BinTree.h"
#include "FileRecordMgr.h"


BinTree::BinTree(ContextIntersect *context) 
:  _context(context), _lastQueryBins(NULL), _haveQueryChr(false), _lastAddBins(NULL),
   _compact(false), _poolUsed(0), _text(false), _textUsed(0), _lastTextOff(0), _keyBlocks(NULL), _stopAtFirst(false) {
}

BinTree::~BinTree() {
	for (mainMapType::iterator it = _mainMap.begin(); it != _mainMap.end(); ++it) {
		delete it->second;
	}
	for (size_t i = 0; i < _poolBlocks.size(); i++) {
		delete [] _poolBlocks[i];
	}
	for (size_t i = 0; i < _textBlocks.size(); i++) {
		delete [] _textBlocks[i];
	}
}

CompactHit *BinTree::nextPoolRecord() {
	if (_poolUsed == _pool.size()) {
		const size_t block = 256;
		CompactHit *b = new CompactHit[block];
		_poolBlocks.push_back(b);
		for (size_t i = 0; i < block; i++) _pool.push_back(b + i);
	}
	return _pool[_poolUsed];
}

// Appends the text the record prints (its coordinates as read) to the arena and returns
// its offset. A zero-length record is printed with its coordinates undone.
int64_t BinTree::appendText(Record *record)
{
	record->undoZeroLength();
	_textScratch.clear();
	record->print(_textScratch);
	record->adjustZeroLength();
	return appendRaw(_textScratch.data(), _textScratch.size());
}

bool BinTree::textEquals(int64_t off, const char *text, size_t len) const
{
	const char *src = _textBlocks[off / TEXT_BLOCK] + off % TEXT_BLOCK;
	uint32_t l;
	memcpy(&l, src, sizeof(l));
	return l == len && memcmp(src + sizeof(l), text, len) == 0;
}

int64_t BinTree::appendRaw(const char *text, size_t len_)
{
	size_t need = sizeof(uint32_t) + len_;
	if (_textBlocks.empty() || _textUsed + need > TEXT_BLOCK) {
		_textBlocks.push_back(new char[need > TEXT_BLOCK ? need : TEXT_BLOCK]);
		_textUsed = 0;
		if (need > TEXT_BLOCK) {
			// an oversized line takes a block of its own
			int64_t off = (int64_t)(_textBlocks.size() - 1) * (int64_t)TEXT_BLOCK;
			_lastTextOff = off;
			uint32_t len = (uint32_t)len_;
			memcpy(_textBlocks.back(), &len, sizeof(len));
			memcpy(_textBlocks.back() + sizeof(len), text, len_);
			_textUsed = TEXT_BLOCK;
			return off;
		}
	}
	int64_t off = (int64_t)(_textBlocks.size() - 1) * (int64_t)TEXT_BLOCK + (int64_t)_textUsed;
	_lastTextOff = off;
	char *dst = _textBlocks.back() + _textUsed;
	uint32_t len = (uint32_t)len_;
	memcpy(dst, &len, sizeof(len));
	memcpy(dst + sizeof(len), text, len_);
	_textUsed += need;
	return off;
}

void BinTree::loadDB()
{
	// Compact when no database record is read beyond its coordinates, file index and
	// zero-length flag, and hits are not sorted (sorting compares whole records).
	_compact = !_context->getSortOutput();
	// -u and -v read only whether a query has any hit (reportOverlapSummary).
	_stopAtFirst = _context->getProgram() == ContextBase::INTERSECT &&
			(_context->getAnyHit() != _context->getNoHit()) &&
			!_context->getWriteCount() && !_context->getWriteCountsPerDatabase() &&
			!_context->getObeySplits();
	for (int i=0; _compact && i < _context->getNumDatabaseFiles(); i++) {
		FileRecordMgr *f = _context->getDatabaseFile(i);
		FileRecordTypeChecker::RECORD_TYPE t = f->getRecordType();
		if (!f->getCoordsOnly() || !(t == FileRecordTypeChecker::BED3_RECORD_TYPE ||
				t == FileRecordTypeChecker::BED4_RECORD_TYPE || t == FileRecordTypeChecker::BED5_RECORD_TYPE ||
				t == FileRecordTypeChecker::BED6_RECORD_TYPE || t == FileRecordTypeChecker::BED_PLUS_RECORD_TYPE ||
				t == FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE)) {
			_compact = false;
		}
	}
	// Text mode: the same, but the database record is printed, so its text is kept.
	_text = false;
	if (!_compact && !_context->getSortOutput() && _context->printOnlyDatabaseFiles()
			&& _context->getNumDatabaseFiles() < 1000) {
		_text = true;
		for (int i=0; _text && i < _context->getNumDatabaseFiles(); i++) {
			FileRecordMgr *f = _context->getDatabaseFile(i);
			FileRecordTypeChecker::RECORD_TYPE t = f->getRecordType();
			if (f->getCoordsOnly() || !(t == FileRecordTypeChecker::BED3_RECORD_TYPE ||
					t == FileRecordTypeChecker::BED4_RECORD_TYPE || t == FileRecordTypeChecker::BED5_RECORD_TYPE ||
					t == FileRecordTypeChecker::BED6_RECORD_TYPE || t == FileRecordTypeChecker::BED_PLUS_RECORD_TYPE ||
					t == FileRecordTypeChecker::BED6_PLUS_RECORD_TYPE)) {
				_text = false;
			}
		}
		if (_text) _compact = true;
	}
	for (int i=0; i < _context->getNumDatabaseFiles(); i++) {
		FileRecordMgr *databaseFile = _context->getDatabaseFile(i);

		Record *record = NULL;
		if (_compact) {
			// Lines of the previous line's chromosome with plain coordinates skip the
			// Record; the first line of a chromosome and every other line take the full path.
			// In text mode a line is also kept as its own text, once the first line of the
			// chromosome has printed to exactly its line text (textLine).
			bool textLine = false;
			string lastChr;
			const string *lastChrPtr = NULL;
			while (!databaseFile->eof()) {
				CHRPOS start = 0, end = 0;
				bool zeroLen = false;
				const char *line = NULL;
				size_t lineLen = 0;
				int got = databaseFile->getNextCoords((_text && !textLine) ? NULL : lastChrPtr, start, end, zeroLen, record,
						&line, &lineLen);
				if (got == 0) continue;
				if (got == 1) {
					int64_t meta = ((int64_t)databaseFile->getFileIdx() << 1) | (zeroLen ? 1 : 0);
					if (_text) {
						meta |= appendRaw(line, lineLen) << 16;
					}
					if (!stageCoords(start, end, meta)) {
						fprintf(stderr, "ERROR: Unable to add record to tree.\n");
						databaseFile->close();
						exit(1);
					}
					continue;
				}
				if (record->isUnmapped()) {
					databaseFile->deleteRecord(record);
					continue;
				}
				_context->testNameConventions(record);
				if (!addRecordToTree(record)) {
					fprintf(stderr, "ERROR: Unable to add record to tree.\n");
					databaseFile->close();
					exit(1);
				}
				if (_text && (lastChrPtr == NULL || lastChr != record->getChrName())) {
					const char *l;
					size_t ll;
					databaseFile->getCurrentLine(l, ll);
					textLine = textEquals(_lastTextOff, l, ll);
				}
				lastChr = record->getChrName();
				lastChrPtr = &lastChr;
				databaseFile->deleteRecord(record);
			}
			continue;
		}
		while (!databaseFile->eof()) {
			record = databaseFile->getNextRecord();
			//In addition to NULL records, we also don't want to add unmapped reads.
			if (record == NULL || record->isUnmapped()) {
				continue;
			}

			_context->testNameConventions(record);

			if (!addRecordToTree(record)) {
				fprintf(stderr, "ERROR: Unable to add record to tree.\n");
				databaseFile->close();
				exit(1);
			}
			// The scratch Record is recycled by the next read.
			if (_compact) {
				databaseFile->deleteRecord(record);
			}
		}
	}
	for (mainMapType::iterator it = _mainMap.begin(); it != _mainMap.end(); ++it) {
		it->second->finalize();
	}
	_haveQueryChr = false;
}

void BinTree::getHits(Record *record, RecordKeyVector &hitSet)
{
	_poolUsed = 0;
	if (record->isUnmapped()) {
		return;
	}
    const ChrRef &chr = record->getChrRef();
	if (!_haveQueryChr || !(chr == _lastQueryChr)) {
		mainMapType::iterator mainIter = _mainMap.find(chr.str());
		_lastQueryChr = chr;
		_haveQueryChr = true;
		_lastQueryBins = (mainIter == _mainMap.end()) ? NULL : mainIter->second;
		for (int l = 0; l < 8; l++) { _curFirst[l] = NULL; _curCnt[l] = 0; _memoBin[l] = -1; }
		_curStart = 0;
	}
	if (_lastQueryBins == NULL) {
		//given chrom not even in map.
		return;
	}

    binNumType startPos = record->getStartPos();
    binNumType endPos = record->getEndPos();

    binNumType startBin = (startPos >> _binFirstShift);
    binNumType endBin = ((endPos-1) >> _binFirstShift);
	if (startPos < _curStart) {
		for (int l = 0; l < 8; l++) { _curFirst[l] = NULL; _curCnt[l] = 0; _memoBin[l] = -1; }
	}
	_curStart = startPos;


	const allBinsType &bins = *_lastQueryBins;

    /* SYNOPSIS:
         1. We loop through each UCSC BIN level for feature A's chrom.
         2. For each BIN, we loop through each B feature and add it to
            hits if it meets all of the user's requests, which include:
               (a) overlap fraction, (b) strandedness, (c) reciprocal overlap
    */
    for (binNumType i = 0; i < NUM_BIN_LEVELS; i++) {
        binNumType offset = _binOffsetsExtended[i];
        for (binNumType j = (startBin+offset); j <= (endBin+offset); j++)  {

        	// move to the next bin if this one is empty
        	const binEntry *iter, *binEnd;
        	bool binSorted;
        	const bool isStartBin = (j == startBin + offset);
        	if (isStartBin) {
        		// the previous query's start bin at this level is often this query's
        		if (_memoBin[i] == j) {
        			if (!_memoFound[i]) continue;
        			iter = _memoFirst[i];
        			binEnd = _memoLast[i];
        			binSorted = _memoSorted[i];
        		} else {
        			bool found = bins.find(j, iter, binEnd, binSorted);
        			_memoBin[i] = j;
        			_memoFound[i] = found;
        			if (!found) continue;
        			_memoFirst[i] = iter;
        			_memoLast[i] = binEnd;
        			_memoSorted[i] = binSorted;
        		}
        	} else if (!bins.find(j, iter, binEnd, binSorted)) {
        		continue;
        	}
        	if (isStartBin) {
        		// skip leading entries already known to end before this query's start
        		const binEntry *first = iter;
        		if (_curFirst[i] == first) iter += _curCnt[i];
        		while (iter != binEnd && iter->end < startPos) iter++;
        		_curFirst[i] = first;
        		_curCnt[i] = (size_t)(iter - first);
        	}

        	for (; iter != binEnd; iter++) {
            	// Record::sameChromIntersects returns false when minEnd < maxStart.
            	if (iter->end < startPos || iter->start > endPos) {
            		// entries after a start past the query's end cannot hit
            		if (binSorted && iter->start > endPos) break;
            		continue;
            	}
            	if (_keyBlocks != NULL && _compact) {
            		bool any = false;
            		for (size_t k = 0; k < _keyBlocks->size(); k++) {
            			if (min((*_keyBlocks)[k].second, iter->end) > max((*_keyBlocks)[k].first, iter->start)) {
            				any = true;
            				break;
            			}
            		}
            		if (!any) continue;
            	}
            	Record *dbRec;
            	if (_compact) {
            		CompactHit *h = nextPoolRecord();
            		if (_text) {
            			int64_t off = iter->meta >> 16;
            			h->setText(_textBlocks[off / TEXT_BLOCK] + off % TEXT_BLOCK);
            			h->setHitFields(_lastQueryBins->chr, iter->start, iter->end,
            					(int)((iter->meta & 0xffff) >> 1), (iter->meta & 1) != 0);
            		} else {
            			h->setHitFields(_lastQueryBins->chr, iter->start, iter->end,
            					(int)(iter->meta >> 1), (iter->meta & 1) != 0);
            		}
            		dbRec = h;
            	} else {
            		dbRec = iter->record;
            	}
            	if (record->intersects(dbRec,
                                       _context->getSameStrand(),
                                       _context->getDiffStrand(),
            			               _context->getOverlapFractionA(),
                                       _context->getOverlapFractionB(),
                                       _context->getReciprocalFraction(),
                                       _context->getEitherFraction(),
                                       _context->getObeySplits()
                                      )
                    )
                {
            		hitSet.push_back(dbRec);
            		if (_compact) _poolUsed++;
            		if (_stopAtFirst) return;
            	}
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
	if (_context->getSortOutput()) {
		hitSet.sortVector();
	}
}

// Stages an entry for the chromosome of the last addRecordToTree call.
bool BinTree::stageCoords(binNumType startPos, binNumType endPos, int64_t meta)
{
	binNumType binNum = getBin(startPos, endPos);
	if (binNum < 0 || binNum >= NUM_BINS) {
		fprintf(stderr, "ERROR: Received illegal bin number %" PRId_BINNUMTYPE " from getBin call.\n"
                                "Maximum values is: %" PRId_BINNUMTYPE "\n"
                                "This typically means that your coordinates are\n"
                                "negative or too large to represent in the data\n"
                                "structure bedtools uses to find intersections.", binNum, NUM_BINS);
		return false;
	}
	_haveQueryChr = false;
	binEntry entry;
	entry.start = startPos;
	entry.end = endPos;
	entry.meta = meta;
	_lastAddBins->stage(binNum, entry);
	return true;
}

bool BinTree::addRecordToTree(Record *record)
{
	// Get chr, bin.
	const string &chr = record->getChrName();
	binNumType startPos = (binNumType)(record->getStartPos());
	binNumType endPos = (binNumType)(record->getEndPos());
	binNumType binNum = getBin(startPos, endPos);

	if (binNum < 0 || binNum >= NUM_BINS) {
		fprintf(stderr, "ERROR: Received illegal bin number %" PRId_BINNUMTYPE " from getBin call.\n"
                                "Maximum values is: %" PRId_BINNUMTYPE "\n"
                                "This typically means that your coordinates are\n"
                                "negative or too large to represent in the data\n"
                                "structure bedtools uses to find intersections.", binNum, NUM_BINS);
		return false;
	}
	if (_lastAddBins == NULL || chr != _lastAddChr) {
		allBinsType *&slot = _mainMap[chr];
		if (slot == NULL) {
			slot = new allBinsType();
			slot->chr = record->getChrRef();
		}
		_lastAddChr = chr;
		_lastAddBins = slot;
	}
	allBinsType *bins = _lastAddBins;
	_haveQueryChr = false;
	binEntry entry;
	entry.start = startPos;
	entry.end = endPos;
	if (_compact) {
		entry.meta = ((int64_t)record->getFileIdx() << 1) | (record->isZeroLength() ? 1 : 0);
		if (_text) {
			entry.meta |= appendText(record) << 16;
		}
	} else {
		entry.record = record;
	}
	bins->stage(binNum, entry);
	return true;
}


BinTree::binNumType BinTree::getBin(const Record *record) const {
	return getBin((binNumType)(record->getStartPos()), (binNumType)(record->getEndPos()));
}

BinTree::binNumType BinTree::getBin(binNumType start, binNumType end) const {
    --end;
    start >>= _binFirstShift;
    end   >>= _binFirstShift;

    for (binNumType i = 0; i < NUM_BIN_LEVELS; ++i) {
        if (start == end) {
        	return _binOffsetsExtended[i] + start;
        }
        start >>= _binNextShift;
        end   >>= _binNextShift;
    }
    //failure
    return -1;
}
