#include "BinTree.h"
#include "FileRecordMgr.h"
#include "Bed3Interval.h"
#include "WorkerPool.h"
#include "ParseTools.h"

#include <climits>
#include <algorithm>
#include <atomic>

// A hit that carries only its B record's coordinates.
class BinTree::CoordRecord final : public Bed3Interval {
public:
	void setCoords(CHRPOS start, CHRPOS end) {
		_startPos = start;
		_endPos = end;
	}
};

// A hit of B loaded as lines, standing in for the line's Record where
// intersect prints or tags hits: it carries the Record's coordinates (as
// adjusted), zero-length state and file index, which output reads, undoes
// and redoes, and prints as the Record would in the same state.
class BinTree::LineHit final : public Bed3Interval {
public:
	explicit LineHit(BinTree *tree) : _tree(tree), _line(NULL) {}
	void set(const BinEntry &entry, int fileIdx) {
		_startPos = entry.start;
		_endPos = entry.end;
		_zeroLength = entry.line->zeroLength;
		_fileIdx = fileIdx;
		_line = entry.line;
	}
	void print(string &outBuf) const override {
		// undone, or never adjusted, a verbatim line's Record prints its text
		if (_line->verbatim && !_zeroLength) {
			outBuf.append(_line->text(), _line->len);
			return;
		}
		// the Record stays adjusted between hits: print it in this state
		Record *rec = _tree->lineRecord(_line);
		bool undo = !_zeroLength && rec->isZeroLength();
		if (undo) {
			rec->undoZeroLength();
		}
		rec->print(outBuf);
		if (undo) {
			rec->adjustZeroLength();
		}
	}
	void print(string &outBuf, CHRPOS start, CHRPOS end) const override {
		_tree->lineRecord(_line)->print(outBuf, start, end);
	}
	void print(string &outBuf, const string &start, const string &end) const override {
		_tree->lineRecord(_line)->print(outBuf, start, end);
	}
private:
	BinTree *_tree;
	LazyLine *_line;
};

// -threads: the parallel load's state (see loadLinesParallel) and findHits'
// per-thread state.
struct BinTree::ChromBins::Piece {
	const BinEntry *recs;
	const uint32_t *keys;
	size_t n;
};

struct BinTree::LoadFormat {
	int numFields;
	char delim;
	bool keepLines; // store each entry's line (!_coordHits)
	bool lineFlags; // set zeroLength and verbatim (_lineHits)
	bool asText;    // the file prints as text (loadLines' asText)
	int fileNum;
};

struct BinTree::LoadBlock {
	// A run of entries on one chromosome, or one line for the main thread.
	struct Run {
		const char *line;      // the run's first line, or the slow line
		size_t lineLen;
		bool slow;
		size_t chrLen;         // the chromosome is line[0, chrLen)
		size_t first;          // entries[first, first + count)
		size_t count;
		unsigned long totalLength;
		bool zeroLength;       // an entry with start == end (stage()'s test)
		binNumType levelLo[NUM_BIN_LEVELS];
		binNumType levelHi[NUM_BIN_LEVELS];
	};
	char *buf;    // bytes whole lines, then a '\n'
	size_t bytes;
	bool firstLineSlow;
	bool parsed;  // guarded by loadLinesParallel's mutex
	vector<BinEntry> entries;
	vector<uint32_t> keys;
	vector<Run> runs;
	vector<char *> arena; // the LazyLines, handed to the tree's arena
	LoadBlock() : buf(NULL), bytes(0), firstLineSlow(false), parsed(false) {}
	~LoadBlock() { delete [] buf; }
};

struct BinTree::ParallelLoad {
	// A piece of a chromosome's records: a run of a block's entries, or
	// (slow != NULL) the entries of slow lines between runs.
	struct SlowPiece {
		vector<BinEntry> recs;
		vector<uint32_t> keys;
	};
	struct PieceRef {
		const BinEntry *recs;
		const uint32_t *keys;
		size_t n;
		SlowPiece *slow;
	};
	map<ChromBins *, vector<PieceRef> > pieces; // each chromosome's, in order
	vector<PieceRef> *pieceList;                // _loadBins'
	vector<LoadBlock *> blocks;                 // holding the runs' entries
	vector<SlowPiece *> slowPieces;
	string slowChr;
	ParallelLoad() : pieceList(NULL) {}
	~ParallelLoad() {
		for (size_t i = 0; i < blocks.size(); i++) {
			delete blocks[i];
		}
		for (size_t i = 0; i < slowPieces.size(); i++) {
			delete slowPieces[i];
		}
	}
};

struct BinTree::QueryState {
	char pad[64]; // apart from the previous thread's
	const ChromBins *rangeBins;
	binNumType rangeStartBin;
	binNumType rangeEndBin;
	int numRanges;
	RecRange ranges[MAX_RANGES];
	QueryState() : rangeBins(NULL), rangeStartBin(0), rangeEndBin(0), numRanges(0) {}
};

// The stand-ins the printing thread serves, written for every query, in a
// block of their own.
struct BinTree::ServeState {
	char pad[64];
	vector<CoordRecord *> coordPool;
	size_t coordUsed;
	vector<LineHit *> linePool;
	size_t lineUsed;
	char pad2[64];
	ServeState() : coordUsed(0), lineUsed(0) {}
	~ServeState() {
		for (size_t i = 0; i < coordPool.size(); i++) {
			delete coordPool[i];
		}
		for (size_t i = 0; i < linePool.size(); i++) {
			delete linePool[i];
		}
	}
};

BinTree::BinTree(ContextIntersect *context) 
:  _context(context),
   _loadBins(NULL),
   _queryBins(NULL),
   _haveQuery(false),
   _rangeBins(NULL),
   _rangeStartBin(0),
   _rangeEndBin(0),
   _numRanges(0),
   _lazy(false),
   _coordHits(false),
   _coordUsed(0),
   _lineHits(false),
   _lineUsed(0),
   _arenaPos(NULL),
   _arenaLeft(0),
   _pool(NULL),
   _parLoad(NULL),
   _serve(NULL) {
}

BinTree::~BinTree() {
	// Records built from lines copied what they needed
	for (size_t i = 0; i < _arenaChunks.size(); i++) {
		delete [] _arenaChunks[i];
	}
	for (size_t i = 0; i < _coordPool.size(); i++) {
		delete _coordPool[i];
	}
	for (size_t i = 0; i < _linePool.size(); i++) {
		delete _linePool[i];
	}
	for (size_t i = 0; i < _queryStates.size(); i++) {
		delete _queryStates[i];
	}
	delete _serve;
}

void BinTree::loadDB()
{
	_lazy = linesLoadable();
	// Of a hit, intersect's output then reads the count (-u, -v, -c) or its
	// start and end (the default report, -wa); -wb, -wo, -wao and -loj
	// print it, -C reads its file, -split its blocks and -sortout sorts by
	// its chromosome. A lazy B record's coordinates are those of its Record
	// (getNextLine adjusts them as coordsValid does), and a stand-in's
	// adjustZeroLength after printing does nothing, since no stored entry
	// has start == end. Such a hit reads only the entry's coordinates, so
	// loadLines keeps no line for it (every entry's line is NULL).
	_coordHits = _lazy && _context->getProgram() == ContextBase::INTERSECT
	             && !_context->getWriteB() && !_context->getWriteOverlap()
	             && !_context->getWriteAllOverlap() && !_context->getLeftJoin()
	             && !_context->getWriteCountsPerDatabase()
	             && !_context->getObeySplits() && !_context->getSortOutput();
	// Where intersect prints or tags hits (-wb, -wo, -wao, -loj, -C), a hit
	// is read for its start, end, zero-length state and file and printed:
	// a line stand-in has all of those.
	_lineHits = _lazy && !_coordHits && _context->getProgram() == ContextBase::INTERSECT
	            && !_context->getObeySplits() && !_context->getSortOutput();
	// with -threads, B text files are parsed on the pool's threads
	if (_lazy && _pool != NULL) {
		bool bulk = true;
		for (int i=0; i < _context->getNumDatabaseFiles(); i++) {
			bulk = bulk && _context->getDatabaseFile(i)->bulkLoadable();
		}
		if (bulk) {
			_parLoad = new ParallelLoad();
		}
	}
	for (int i=0; i < _context->getNumDatabaseFiles(); i++) {
		FileRecordMgr *databaseFile = _context->getDatabaseFile(i);
		if (_lazy) {
			if (_parLoad != NULL) {
				loadLinesParallel(i, databaseFile);
			} else {
				loadLines(i, databaseFile);
			}
			continue;
		}

		Record *record = NULL;
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
		}
	}
	if (_parLoad != NULL) {
		finalizeParallel();
	}
	// nothing more is staged: free each chromosome's blocks as it is indexed,
	// so that the memory serves the next chromosome's recs
	_stagePool.release();
	for (mainMapType::iterator it = _mainMap.begin(); it != _mainMap.end(); ++it) {
		it->second.finalize(_stagePool);
		_stagePool.release();
	}
	_haveQuery = false;
	_rangeBins = NULL;
}

// B can be loaded as lines when every database file allows it and getHits
// decides overlap from the stored coordinates alone, so that only hits need
// their Records: no strand test and no overlap fraction. The query's
// chromosome id is -1 too, since a genome file stops the files allowing it;
// getHits builds each candidate's Record should that ever not hold, except
// for coordinate stand-ins (_coordHits), which keep no line to build from.
bool BinTree::linesLoadable()
{
	if (_context->getSameStrand() || _context->getDiffStrand()
	    || _context->getOverlapFractionA() != 0.0 || _context->getOverlapFractionB() != 0.0) {
		return false;
	}
	// LazyLine::fileNum has 30 bits
	if (_context->getNumDatabaseFiles() >= (1 << 30)) {
		return false;
	}
	for (int i=0; i < _context->getNumDatabaseFiles(); i++) {
		if (!_context->getDatabaseFile(i)->linesLoadable()) {
			return false;
		}
	}
	return true;
}

// The eager loop of loadDB, keeping each record's line in place of its
// Record; every check and message comes at the same line in the same order.
void BinTree::loadLines(int fileNum, FileRecordMgr *databaseFile)
{
	string chr;
	CHRPOS start = 0;
	CHRPOS end = 0;
	unsigned flags = 0;
	bool asText = false;
	if (_lineHits) {
		if ((int)_lineFileIdx.size() <= fileNum) {
			_lineFileIdx.resize(fileNum + 1, -1);
		}
		_lineFileIdx[fileNum] = databaseFile->lineFileIdx();
		asText = databaseFile->linesPrintAsText();
	}
	while (!databaseFile->eof()) {
		if (!databaseFile->getNextLine(chr, start, end, _lineHits ? &flags : NULL)) {
			continue;
		}
		LazyLine *line = NULL;
		if (!_coordHits) {
			line = storeLine(databaseFile->currentLine(), databaseFile->currentLineLength(), fileNum,
			                 (flags & FileRecordMgr::LINE_ZERO_LENGTH) != 0,
			                 asText && (flags & FileRecordMgr::LINE_CANONICAL) != 0);
		}
		if (!_context->testNameConventionsQuietly(chr)) {
			if (line != NULL) {
				_context->testNameConventions(lineRecord(line));
			} else {
				_context->testNameConventions(databaseFile->recordFromLine(databaseFile->currentLine(),
				                                                           databaseFile->currentLineLength()));
			}
		}
		BinEntry entry;
		entry.start = start;
		entry.end = end;
		entry.line = line;
		if (!addEntryToTree(chr, entry)) {
			fprintf(stderr, "ERROR: Unable to add record to tree.\n");
			databaseFile->close();
			exit(1);
		}
	}
}

BinTree::LazyLine *BinTree::storeLine(const char *text, size_t len, int fileNum, bool zeroLength, bool verbatim)
{
	bool keepText = (len <= 0xffffffffu);
	size_t need = sizeof(LazyLine) + (keepText ? (len + 7) & ~(size_t)7 : 0);
	if (need > _arenaLeft) {
		size_t size = need > ARENA_CHUNK ? need : ARENA_CHUNK;
		_arenaPos = new char[size];
		_arenaChunks.push_back(_arenaPos);
		_arenaLeft = size;
	}
	LazyLine *line = reinterpret_cast<LazyLine *>(_arenaPos);
	_arenaPos += need;
	_arenaLeft -= need;
	line->rec = NULL;
	line->len = keepText ? (uint32_t)len : 0;
	line->fileNum = (uint32_t)fileNum;
	line->zeroLength = zeroLength;
	line->verbatim = keepText && verbatim;
	if (keepText) {
		memcpy(line + 1, text, len);
	} else {
		line->rec = _context->getDatabaseFile(fileNum)->recordFromLine(text, len);
	}
	return line;
}

Record *BinTree::buildLineRecord(LazyLine *line)
{
	line->rec = _context->getDatabaseFile(line->fileNum)->recordFromLine(line->text(), line->len);
	return line->rec;
}

void BinTree::getHits(Record *record, RecordKeyVector &hitSet)
{
	// the previous query's stand-ins are free once its hits are cleared
	_coordUsed = 0;
	_lineUsed = 0;
	if (record->isUnmapped()) {
		return;
	}
    const string &chr = record->getChrName();
	if (!_haveQuery || chr != _queryChr) {
		mainMapType::const_iterator mainIter = _mainMap.find(chr);
		_queryBins = (mainIter == _mainMap.end()) ? NULL : &mainIter->second;
		_queryChr = chr;
		_haveQuery = true;
	}
	if (_queryBins == NULL) {
		//given chrom not even in map.
		return;
	}

    binNumType startPos = record->getStartPos();
    binNumType endPos = record->getEndPos();

    binNumType startBin = (startPos >> _binFirstShift);
    binNumType endBin = ((endPos-1) >> _binFirstShift);


	const ChromBins &bins = *_queryBins;

    /* SYNOPSIS:
         1. We loop through each UCSC BIN level for feature A's chrom.
         2. For each BIN, we loop through each B feature and add it to
            hits if it meets all of the user's requests, which include:
               (a) overlap fraction, (b) strandedness, (c) reciprocal overlap
       Bins are visited by global bin number, ascending within a level; only
       the bins inside a segment can hold records, and a segment's bins from
       lo to hi hold one run of recs. The runs depend on the chromosome and
       the level-0 bins alone, so they are found again only when those change.
    */
	if (_rangeBins != _queryBins || startBin != _rangeStartBin || endBin != _rangeEndBin) {
		const BinSegment *segBegin = bins.segments.data();
		const BinSegment *segEnd = segBegin + bins.segments.size();
		const size_t *binStarts = bins.binStarts.data();
		const BinEntry *recs = bins.recs;
		binNumType levelStart = startBin;
		binNumType levelEnd = endBin;
		int numRanges = 0;
		for (binNumType i = 0; i < NUM_BIN_LEVELS; i++) {
			binNumType offset = _binOffsetsExtended[i];
			binNumType first = levelStart + offset;
			binNumType last = levelEnd + offset;
			for (const BinSegment *seg = segBegin; seg != segEnd; ++seg) {
				if (seg->hi < first) {
					continue;
				}
				if (seg->lo > last) {
					break;
				}
				binNumType lo = first > seg->lo ? first : seg->lo;
				binNumType hi = last < seg->hi ? last : seg->hi;
				const size_t *slot = binStarts + seg->base + (size_t)(lo - seg->lo);
				if (slot[0] < slot[hi - lo + 1]) {
					_ranges[numRanges].begin = recs + slot[0];
					_ranges[numRanges].end = recs + slot[hi - lo + 1];
					numRanges++;
				}
			}
			levelStart >>= _binNextShift;
			levelEnd >>= _binNextShift;
		}
		_numRanges = numRanges;
		_rangeBins = _queryBins;
		_rangeStartBin = startBin;
		_rangeEndBin = endBin;
	}
	const RecRange *rangeEnd = _ranges + _numRanges;

	// Without strand or fraction requirements, Record::intersects reduces to
	// the coordinate test below. It is decided from the stored coordinates when
	// they are the records' own (coordsStable) and the chromosome test is a
	// name compare (no chromosome id on the query), which the map lookup has
	// already passed; B's unmapped records were never loaded and an unmapped
	// query returned above. For coordinate stand-ins (whose entries keep no
	// line) the test is exact for any query: their B file has no genome file,
	// so its Records would have chromosome id -1 and sameChrom would compare
	// names, and getNextLine adjusts start == end, so coordsStable holds.
	bool sameStrand = _context->getSameStrand();
	bool diffStrand = _context->getDiffStrand();
	float overlapFractionA = _context->getOverlapFractionA();
	float overlapFractionB = _context->getOverlapFractionB();
	bool coordsOnly = _coordHits
	                  || (bins.coordsStable && !sameStrand && !diffStrand
	                      && overlapFractionA == 0.0 && overlapFractionB == 0.0
	                      && record->getChromId() == -1);
	bool localZeroLen = (endPos - startPos == 0);

	if (coordsOnly) {
		// Record::sameChromIntersects with both fractions zero
		for (const RecRange *range = _ranges; range != rangeEnd; ++range) {
			for (const BinEntry *dbEntry = range->begin, *runEnd = range->end; dbEntry != runEnd; ++dbEntry) {
				CHRPOS maxStart = max(startPos, dbEntry->start);
				CHRPOS minEnd = min(endPos, dbEntry->end);
				if (minEnd < maxStart) {
					continue;
				}
				if (minEnd == maxStart && !localZeroLen && dbEntry->end - dbEntry->start != 0) {
					continue;
				}
				if (_coordHits) {
					CoordRecord *hit;
					if (_coordUsed < _coordPool.size()) {
						hit = _coordPool[_coordUsed];
					} else {
						hit = new CoordRecord();
						_coordPool.push_back(hit);
					}
					_coordUsed++;
					hit->setCoords(dbEntry->start, dbEntry->end);
					hitSet.push_back(hit);
					continue;
				}
				if (_lineHits) {
					LineHit *hit;
					if (_lineUsed < _linePool.size()) {
						hit = _linePool[_lineUsed];
					} else {
						hit = new LineHit(this);
						_linePool.push_back(hit);
					}
					_lineUsed++;
					hit->set(*dbEntry, _lineFileIdx[dbEntry->line->fileNum]);
					hitSet.push_back(hit);
					continue;
				}
				hitSet.push_back(_lazy ? lineRecord(dbEntry->line) : dbEntry->rec);
			}
		}
	} else {
		for (const RecRange *range = _ranges; range != rangeEnd; ++range) {
			for (const BinEntry *dbEntry = range->begin, *runEnd = range->end; dbEntry != runEnd; ++dbEntry) {
				Record *dbRec = _lazy ? lineRecord(dbEntry->line) : dbEntry->rec;
				if (record->intersects(dbRec,
				                       sameStrand,
				                       diffStrand,
				                       overlapFractionA,
				                       overlapFractionB,
				                       _context->getReciprocalFraction(),
				                       _context->getEitherFraction(),
				                       _context->getObeySplits()
				                      )
				    )
				{
					hitSet.push_back(dbRec);
				}
			}
		}
	}
	if (_context->getSortOutput()) {
		hitSet.sortVector();
	}
}

bool BinTree::addRecordToTree(Record *record)
{
	BinEntry entry;
	entry.start = record->getStartPos();
	entry.end = record->getEndPos();
	entry.rec = record;
	return addEntryToTree(record->getChrName(), entry);
}

bool BinTree::addEntryToTree(const string &chr, const BinEntry &entry)
{
	// Get bin.
	binNumType startPos = (binNumType)(entry.start);
	binNumType endPos = (binNumType)(entry.end);
	int level = -1;
	binNumType binNum = getBinLevel(startPos, endPos, level);

	if (binNum < 0 || binNum >= NUM_BINS) {
		fprintf(stderr, "ERROR: Received illegal bin number %" PRId_BINNUMTYPE " from getBin call.\n"
                                "Maximum values is: %" PRId_BINNUMTYPE "\n"
                                "This typically means that your coordinates are\n"
                                "negative or too large to represent in the data\n"
                                "structure bedtools uses to find intersections.", binNum, NUM_BINS);
		return false;
	}
	if (_loadBins == NULL || chr != _loadChr) {
		// The first run of a chromosome's records has ended: index them now,
		// so that the next chromosome stages into the same blocks. Later runs
		// wait for loadDB, so that interleaved chromosomes cost one more
		// pass at most.
		if (_loadBins != NULL && _loadBins->recs == NULL) {
			_loadBins->finalize(_stagePool);
		}
		_loadBins = &_mainMap[chr];
		_loadChr = chr;
	}
	_loadBins->stage(entry, binNum, level, _stagePool);
	return true;
}


BinTree::binNumType BinTree::getBin(const Record *record) const {
	return getBin((binNumType)(record->getStartPos()), (binNumType)(record->getEndPos()));
}

BinTree::binNumType BinTree::getBin(binNumType start, binNumType end) const {
	int level;
	return getBinLevel(start, end, level);
}

// The global bin number of [start, end) and the level it lies on (-1 and -1
// when no level holds it).
BinTree::binNumType BinTree::getBinLevel(binNumType start, binNumType end, int &level) {
    --end;
    start >>= _binFirstShift;
    end   >>= _binFirstShift;

    for (binNumType i = 0; i < NUM_BIN_LEVELS; ++i) {
        if (start == end) {
        	level = (int)i;
        	return _binOffsetsExtended[i] + start;
        }
        start >>= _binNextShift;
        end   >>= _binNextShift;
    }
    //failure
    level = -1;
    return -1;
}

BinTree::StageBlock BinTree::StagePool::take(int sizeClass) {
	StageBlock block;
	if (!blocks[sizeClass].empty()) {
		block = blocks[sizeClass].back();
		blocks[sizeClass].pop_back();
		return block;
	}
	size_t size = STAGE_BLOCK_MIN << sizeClass;
	block.recs = new BinEntry[size];
	block.keys = new uint32_t[size];
	block.sizeClass = sizeClass;
	return block;
}

void BinTree::StagePool::release() {
	for (int c = 0; c < STAGE_SIZE_CLASSES; c++) {
		for (size_t b = 0; b < blocks[c].size(); b++) {
			delete [] blocks[c][b].recs;
			delete [] blocks[c][b].keys;
		}
		vector<StageBlock>().swap(blocks[c]);
	}
}

BinTree::ChromBins::ChromBins()
:  stageRecPos(NULL),
   stageKeyPos(NULL),
   stageLeft(0),
   numStaged(0),
   recs(NULL),
   numRecs(0),
   coordsStable(true) {
	for (int i = 0; i < NUM_BIN_LEVELS; i++) {
		levelLo[i] = 0;
		levelHi[i] = -1;
	}
}

BinTree::ChromBins::~ChromBins() {
	for (size_t b = 0; b < stageBlocks.size(); b++) {
		delete [] stageBlocks[b].recs;
		delete [] stageBlocks[b].keys;
	}
	delete [] recs;
}

void BinTree::ChromBins::addStageBlock(StagePool &pool) {
	int sizeClass = 0;
	if (!stageBlocks.empty()) {
		sizeClass = stageBlocks.back().sizeClass + 1;
		if (sizeClass >= STAGE_SIZE_CLASSES) {
			sizeClass = STAGE_SIZE_CLASSES - 1;
		}
	}
	StageBlock block = pool.take(sizeClass);
	stageBlocks.push_back(block);
	stageRecPos = block.recs;
	stageKeyPos = block.keys;
	stageLeft = STAGE_BLOCK_MIN << sizeClass;
}

// binNum is in [0, NUM_BINS), checked by addRecordToTree.
void BinTree::ChromBins::stage(const BinEntry &entry, binNumType binNum, int level, StagePool &pool) {
	if (levelHi[level] < levelLo[level]) {
		levelLo[level] = binNum;
		levelHi[level] = binNum;
	} else if (binNum < levelLo[level]) {
		levelLo[level] = binNum;
	} else if (binNum > levelHi[level]) {
		levelHi[level] = binNum;
	}
	if (stageLeft == 0) {
		addStageBlock(pool);
	}
	*stageRecPos++ = entry;
	if (entry.start == entry.end) {
		coordsStable = false;
	}
	*stageKeyPos++ = (uint32_t)binNum | ((uint32_t)level << STAGE_LEVEL_SHIFT);
	stageLeft--;
	numStaged++;
}

void BinTree::ChromBins::finalize(StagePool &pool) {
	if (numStaged == 0) {
		return;
	}
	// Each level's range of bins; ranges that overlap or touch are merged, so
	// every level's range lies inside exactly one segment.
	vector<BinSegment> ranges;
	for (int i = 0; i < NUM_BIN_LEVELS; i++) {
		if (levelHi[i] >= levelLo[i]) {
			BinSegment r;
			r.lo = levelLo[i];
			r.hi = levelHi[i];
			r.base = 0;
			ranges.push_back(r);
		}
	}
	// insertion sort by lo: at most NUM_BIN_LEVELS ranges
	for (size_t a = 1; a < ranges.size(); a++) {
		BinSegment r = ranges[a];
		size_t b = a;
		while (b > 0 && ranges[b - 1].lo > r.lo) {
			ranges[b] = ranges[b - 1];
			b--;
		}
		ranges[b] = r;
	}
	segments.clear();
	for (size_t a = 0; a < ranges.size(); a++) {
		if (!segments.empty() && ranges[a].lo <= segments.back().hi + 1) {
			if (ranges[a].hi > segments.back().hi) {
				segments.back().hi = ranges[a].hi;
			}
		} else {
			segments.push_back(ranges[a]);
		}
	}
	size_t numBins = 0;
	for (size_t a = 0; a < segments.size(); a++) {
		segments[a].base = numBins;
		numBins += (size_t)(segments[a].hi - segments[a].lo + 1);
	}

	// the slot of each level's lowest bin, less that bin, so a staged record's
	// slot is levelSlot[level] + bin
	int64_t levelSlot[NUM_BIN_LEVELS];
	for (int i = 0; i < NUM_BIN_LEVELS; i++) {
		levelSlot[i] = 0;
		if (levelHi[i] < levelLo[i]) {
			continue;
		}
		for (size_t a = 0; a < segments.size(); a++) {
			if (segments[a].lo <= levelLo[i] && levelHi[i] <= segments[a].hi) {
				levelSlot[i] = (int64_t)segments[a].base - segments[a].lo;
				break;
			}
		}
	}

	// stable counting sort by slot of the records already in recs, in their
	// order, then of the staged ones, block by block in insertion order;
	// each bin's records in recs are in insertion order and were all
	// inserted before the staged ones
	const uint32_t binMask = (1u << STAGE_LEVEL_SHIFT) - 1;
	size_t numBlocks = stageBlocks.size();
	binStarts.assign(numBins + 1, 0);
	for (size_t r = 0; r < numRecs; r++) {
		int level = -1;
		binNumType bin = getBinLevel((binNumType)recs[r].start, (binNumType)recs[r].end, level);
		binStarts[(size_t)(levelSlot[level] + bin) + 1]++;
	}
	for (size_t b = 0; b < numBlocks; b++) {
		const uint32_t *keys = stageBlocks[b].keys;
		size_t n = (STAGE_BLOCK_MIN << stageBlocks[b].sizeClass) - (b + 1 == numBlocks ? stageLeft : 0);
		for (size_t r = 0; r < n; r++) {
			uint32_t key = keys[r];
			size_t slot = (size_t)(levelSlot[key >> STAGE_LEVEL_SHIFT] + (int64_t)(key & binMask));
			binStarts[slot + 1]++;
		}
	}
	for (size_t k = 0; k < numBins; k++) {
		binStarts[k + 1] += binStarts[k];
	}
	vector<size_t> cursor(binStarts.begin(), binStarts.end() - 1);
	size_t total = numRecs + numStaged;
	BinEntry *sorted = new BinEntry[total];
	for (size_t r = 0; r < numRecs; r++) {
		int level = -1;
		binNumType bin = getBinLevel((binNumType)recs[r].start, (binNumType)recs[r].end, level);
		sorted[cursor[(size_t)(levelSlot[level] + bin)]++] = recs[r];
	}
	for (size_t b = 0; b < numBlocks; b++) {
		const uint32_t *keys = stageBlocks[b].keys;
		const BinEntry *staged = stageBlocks[b].recs;
		size_t n = (STAGE_BLOCK_MIN << stageBlocks[b].sizeClass) - (b + 1 == numBlocks ? stageLeft : 0);
		for (size_t r = 0; r < n; r++) {
			uint32_t key = keys[r];
			size_t slot = (size_t)(levelSlot[key >> STAGE_LEVEL_SHIFT] + (int64_t)(key & binMask));
			sorted[cursor[slot]++] = staged[r];
		}
	}
	delete [] recs;
	recs = sorted;
	numRecs = total;

	for (size_t b = 0; b < numBlocks; b++) {
		pool.give(stageBlocks[b]);
	}
	stageBlocks.clear();
	stageRecPos = NULL;
	stageKeyPos = NULL;
	stageLeft = 0;
	numStaged = 0;
}

// ---------------------------------------------------------------------------
// -threads: loading B on several threads (loadLinesParallel)
//
// The main thread reads the file's bytes in blocks of whole lines, which the
// pool's threads parse; the main thread then takes each block's lines in
// order. A worker turns each line that getNextLine would read as an entry,
// with no message, and that addEntryToTree takes into an entry and its key
// (and its LazyLine, in the block's own arena, unless _coordHits), runs of
// consecutive such lines on one chromosome becoming one piece of that
// chromosome's records; the main thread does per run what loadLines does per
// line that a run can see (the line count, the name-convention test, which a
// run's first line decides, the record length sum). Every other line (blank,
// a header, a wrong field count, a NUL, coordinates that are not 1 to 18
// digits, end before start, an illegal bin) goes to the main thread's copy of
// loadLines (loadSlowLine) at its turn, so every message and error comes at
// the same line as on one thread. Each chromosome is indexed once at the end
// (finalizePieces), from its pieces in insertion order, which leaves each bin
// in insertion order as finalize does; chromosomes are indexed on the threads.

void BinTree::setPool(WorkerPool *pool)
{
	_pool = pool;
	for (int t = (int)_queryStates.size(); t < pool->size(); t++) {
		_queryStates.push_back(new QueryState());
	}
	if (_serve == NULL) {
		_serve = new ServeState();
	}
}

void BinTree::loadLinesParallel(int fileNum, FileRecordMgr *databaseFile)
{
	static const size_t BLOCK_BYTES = 1 << 20;
	LoadFormat format;
	format.numFields = databaseFile->getNumFields();
	format.delim = databaseFile->textReader()->getDelimChar();
	format.keepLines = !_coordHits;
	format.lineFlags = _lineHits;
	format.asText = false;
	format.fileNum = fileNum;
	if (_lineHits) {
		if ((int)_lineFileIdx.size() <= fileNum) {
			_lineFileIdx.resize(fileNum + 1, -1);
		}
		_lineFileIdx[fileNum] = databaseFile->lineFileIdx();
		format.asText = databaseFile->linesPrintAsText();
	}
	const size_t maxAhead = 4 * (size_t)_pool->size() + 4;

	// blocks[nextParse..] wait for a thread; the main thread alone appends
	// to blocks and takes them in order from nextConsume
	vector<LoadBlock *> blocks;
	size_t nextParse = 0;
	size_t nextConsume = 0;
	bool done = false;
	std::mutex mutex;
	std::condition_variable work;
	std::condition_variable parsed;
	_pool->start([&](int) {
		for (;;) {
			LoadBlock *block;
			{
				std::unique_lock<std::mutex> lock(mutex);
				while (nextParse == blocks.size() && !done) {
					work.wait(lock);
				}
				if (nextParse == blocks.size()) {
					return;
				}
				block = blocks[nextParse++];
			}
			parseBlock(*block, format);
			{
				std::lock_guard<std::mutex> lock(mutex);
				block->parsed = true;
			}
			parsed.notify_all();
		}
	});

	string carry; // the start of a line the last block cut off
	bool atEnd = false;
	bool firstBlock = true;
	for (;;) {
		if (!atEnd && blocks.size() - nextConsume < maxAhead) {
			// the next block: whole lines, all that is left at the end
			size_t cap = max(BLOCK_BYTES, 2 * carry.size());
			char *buf = new char[cap + 1];
			size_t size = carry.size();
			memcpy(buf, carry.data(), size);
			size_t cut = 0;
			for (;;) {
				while (size < cap && !atEnd) {
					size_t n = databaseFile->readRaw(buf + size, cap - size);
					atEnd = (n == 0);
					size += n;
				}
				if (atEnd) {
					cut = size;
					break;
				}
				cut = size;
				while (cut > 0 && buf[cut - 1] != '\n') {
					cut--;
				}
				if (cut > 0) {
					break;
				}
				// a line longer than the block
				char *bigger = new char[2 * cap + 1];
				memcpy(bigger, buf, size);
				delete [] buf;
				buf = bigger;
				cap *= 2;
			}
			carry.assign(buf + cut, size - cut);
			if (cut == 0) {
				delete [] buf;
				continue;
			}
			buf[cut] = '\n';
			LoadBlock *block = new LoadBlock();
			block->buf = buf;
			block->bytes = cut;
			block->firstLineSlow = firstBlock && databaseFile->textReader()->firstLineIsHeader();
			firstBlock = false;
			{
				std::lock_guard<std::mutex> lock(mutex);
				blocks.push_back(block);
			}
			work.notify_one();
			continue;
		}
		if (nextConsume == blocks.size()) {
			break;
		}
		LoadBlock *block = blocks[nextConsume];
		{
			std::unique_lock<std::mutex> lock(mutex);
			if (!block->parsed && nextParse < blocks.size()) {
				// help: parse the next block waiting
				LoadBlock *mine = blocks[nextParse++];
				lock.unlock();
				parseBlock(*mine, format);
				lock.lock();
				mine->parsed = true;
				continue;
			}
			while (!block->parsed) {
				parsed.wait(lock);
			}
		}
		consumeBlock(*block, fileNum, databaseFile);
		nextConsume++;
	}
	{
		std::lock_guard<std::mutex> lock(mutex);
		done = true;
	}
	work.notify_all();
	_pool->wait();
}

// On any thread: the block's lines as loadLinesParallel describes; it
// writes only the block, and reads the block, the format and the tree's
// constants.
void BinTree::parseBlock(LoadBlock &block, const LoadFormat &format) const
{
	const int numFields = format.numFields;
	const char delim = format.delim;
	vector<int> pos(numFields + 1);
	char *arenaPos = NULL;
	size_t arenaLeft = 0;
	const size_t noRun = (size_t)-1;
	size_t runIdx = noRun;
	char *p = block.buf;
	char *blockEnd = block.buf + block.bytes;
	bool firstLine = block.firstLineSlow;
	block.entries.reserve(block.bytes / 24 + 16);
	block.keys.reserve(block.bytes / 24 + 16);
	while (p < blockEnd) {
		// the line as getLine returns it: up to the newline, '\r's before
		// it stripped, a NUL after them
		char *nl = static_cast<char *>(memchr(p, '\n', (size_t)(blockEnd - p) + 1));
		char *next = nl + 1;
		size_t n = (size_t)(nl - p);
		while (n > 0 && p[n - 1] == '\r') {
			n--;
		}
		p[n] = '\0';

		// readEntry's tests: not blank, not a header, the file's field
		// count before any NUL (a line with a NUL goes to the main thread)
		bool clean = !firstLine && n > 0 && n < (size_t)INT_MAX && !isHeaderLine(p, n);
		firstLine = false;
		if (clean) {
			int nf = 1;
			for (size_t i = 0; i < n; i++) {
				char c = p[i];
				if (c == delim) {
					if (nf <= numFields) {
						pos[nf] = (int)i;
					}
					nf++;
				} else if (c == '\0') {
					nf = -1;
					break;
				}
			}
			clean = (nf == numFields);
		}
		// getNextLine's coordinates, without its slow path and its error
		CHRPOS start = 0;
		CHRPOS end = 0;
		bool canonical = false;
		if (clean) {
			pos[0] = -1;
			pos[numFields] = (int)n;
			const char *s1 = p + pos[1] + 1;
			size_t n1 = (size_t)(pos[2] - pos[1] - 1);
			const char *s2 = p + pos[2] + 1;
			size_t n2 = (size_t)(pos[3] - pos[2] - 1);
			clean = str2chrPosDigits(s1, n1, start) && str2chrPosDigits(s2, n2, end) && end >= start;
			canonical = clean && (s1[0] != '0' || n1 == 1) && (s2[0] != '0' || n2 == 1);
		}
		bool zeroLength = (start == end);
		if (clean && zeroLength) {
			start--;
			end++;
		}
		// addEntryToTree's bin
		int level = -1;
		binNumType binNum = -1;
		if (clean) {
			binNum = getBinLevel((binNumType)start, (binNumType)end, level);
			clean = (binNum >= 0 && binNum < NUM_BINS);
		}
		if (!clean) {
			LoadBlock::Run run;
			run.line = p;
			run.lineLen = n;
			run.slow = true;
			run.chrLen = 0;
			run.first = block.entries.size();
			run.count = 0;
			run.totalLength = 0;
			run.zeroLength = false;
			block.runs.push_back(run);
			runIdx = noRun;
			p = next;
			continue;
		}

		size_t chrLen = (size_t)pos[1];
		if (runIdx == noRun || block.runs[runIdx].chrLen != chrLen
		    || memcmp(block.runs[runIdx].line, p, chrLen) != 0) {
			LoadBlock::Run run;
			run.line = p;
			run.lineLen = n;
			run.slow = false;
			run.chrLen = chrLen;
			run.first = block.entries.size();
			run.count = 0;
			run.totalLength = 0;
			run.zeroLength = false;
			for (int i = 0; i < NUM_BIN_LEVELS; i++) {
				run.levelLo[i] = 0;
				run.levelHi[i] = -1;
			}
			block.runs.push_back(run);
			runIdx = block.runs.size() - 1;
		}
		LoadBlock::Run &run = block.runs[runIdx];

		BinEntry entry;
		entry.start = start;
		entry.end = end;
		entry.line = NULL;
		if (format.keepLines) {
			// storeLine, into the block's arena
			size_t need = sizeof(LazyLine) + ((n + 7) & ~(size_t)7);
			if (need > arenaLeft) {
				size_t size = max(need, 2 * (size_t)(blockEnd - p) + 4096);
				arenaPos = new char[size];
				block.arena.push_back(arenaPos);
				arenaLeft = size;
			}
			LazyLine *line = reinterpret_cast<LazyLine *>(arenaPos);
			arenaPos += need;
			arenaLeft -= need;
			line->rec = NULL;
			line->len = (uint32_t)n;
			line->fileNum = (uint32_t)format.fileNum;
			line->zeroLength = format.lineFlags && zeroLength;
			line->verbatim = format.lineFlags && format.asText && canonical;
			memcpy(line + 1, p, n);
			entry.line = line;
		}
		block.entries.push_back(entry);
		block.keys.push_back((uint32_t)binNum | ((uint32_t)level << STAGE_LEVEL_SHIFT));
		run.count++;
		run.totalLength += (unsigned long)(end - start);
		if (start == end) {
			run.zeroLength = true;
		}
		if (run.levelHi[level] < run.levelLo[level]) {
			run.levelLo[level] = binNum;
			run.levelHi[level] = binNum;
		} else if (binNum < run.levelLo[level]) {
			run.levelLo[level] = binNum;
		} else if (binNum > run.levelHi[level]) {
			run.levelHi[level] = binNum;
		}
		p = next;
	}
}

void BinTree::setLoadChrom(const char *chr, size_t len)
{
	if (_loadBins == NULL || _parLoad->pieceList == NULL || _loadChr.size() != len
	    || memcmp(_loadChr.data(), chr, len) != 0) {
		_loadChr.assign(chr, len);
		_loadBins = &_mainMap[_loadChr];
		_parLoad->pieceList = &_parLoad->pieces[_loadBins];
	}
}

// On the main thread, in order: what loadLines does for the block's lines.
void BinTree::consumeBlock(LoadBlock &block, int fileNum, FileRecordMgr *databaseFile)
{
	ParallelLoad &load = *_parLoad;
	bool asText = _lineHits && databaseFile->linesPrintAsText();
	for (size_t i = 0; i < block.runs.size(); i++) {
		const LoadBlock::Run &run = block.runs[i];
		if (run.slow) {
			loadSlowLine(run.line, run.lineLen, fileNum, databaseFile, asText);
			continue;
		}
		setLoadChrom(run.line, run.chrLen);
		databaseFile->skipLines(run.count, run.totalLength);
		// the run's later lines pass, as its first leaves the conventions
		// set or the warning given
		if (!_context->testNameConventionsQuietly(_loadChr)) {
			if (block.entries[run.first].line != NULL) {
				_context->testNameConventions(lineRecord(block.entries[run.first].line));
			} else {
				_context->testNameConventions(databaseFile->recordFromLine(run.line, run.lineLen));
			}
		}
		ChromBins &bins = *_loadBins;
		for (int l = 0; l < NUM_BIN_LEVELS; l++) {
			if (run.levelHi[l] < run.levelLo[l]) {
				continue;
			}
			if (bins.levelHi[l] < bins.levelLo[l]) {
				bins.levelLo[l] = run.levelLo[l];
				bins.levelHi[l] = run.levelHi[l];
			} else {
				bins.levelLo[l] = min(bins.levelLo[l], run.levelLo[l]);
				bins.levelHi[l] = max(bins.levelHi[l], run.levelHi[l]);
			}
		}
		if (run.zeroLength) {
			bins.coordsStable = false;
		}
		ParallelLoad::PieceRef piece;
		piece.recs = block.entries.data() + run.first;
		piece.keys = block.keys.data() + run.first;
		piece.n = run.count;
		piece.slow = NULL;
		load.pieceList->push_back(piece);
	}
	_arenaChunks.insert(_arenaChunks.end(), block.arena.begin(), block.arena.end());
	block.arena.clear();
	delete [] block.buf;
	block.buf = NULL;
	load.blocks.push_back(&block);
}

// loadLines and addEntryToTree for one line a worker left to the main thread.
void BinTree::loadSlowLine(const char *text, size_t len, int fileNum, FileRecordMgr *databaseFile, bool asText)
{
	ParallelLoad &load = *_parLoad;
	string &chr = load.slowChr;
	CHRPOS start = 0;
	CHRPOS end = 0;
	unsigned flags = 0;
	if (!databaseFile->takeNextLine(text, len, chr, start, end, _lineHits ? &flags : NULL)) {
		return;
	}
	LazyLine *line = NULL;
	if (!_coordHits) {
		line = storeLine(databaseFile->currentLine(), databaseFile->currentLineLength(), fileNum,
		                 (flags & FileRecordMgr::LINE_ZERO_LENGTH) != 0,
		                 asText && (flags & FileRecordMgr::LINE_CANONICAL) != 0);
	}
	if (!_context->testNameConventionsQuietly(chr)) {
		if (line != NULL) {
			_context->testNameConventions(lineRecord(line));
		} else {
			_context->testNameConventions(databaseFile->recordFromLine(databaseFile->currentLine(),
			                                                           databaseFile->currentLineLength()));
		}
	}
	int level = -1;
	binNumType binNum = getBinLevel((binNumType)start, (binNumType)end, level);
	if (binNum < 0 || binNum >= NUM_BINS) {
		fprintf(stderr, "ERROR: Received illegal bin number %" PRId_BINNUMTYPE " from getBin call.\n"
                                "Maximum values is: %" PRId_BINNUMTYPE "\n"
                                "This typically means that your coordinates are\n"
                                "negative or too large to represent in the data\n"
                                "structure bedtools uses to find intersections.", binNum, NUM_BINS);
		fprintf(stderr, "ERROR: Unable to add record to tree.\n");
		databaseFile->close();
		exit(1);
	}
	setLoadChrom(chr.data(), chr.size());
	vector<ParallelLoad::PieceRef> &list = *load.pieceList;
	if (list.empty() || list.back().slow == NULL) {
		ParallelLoad::PieceRef piece;
		piece.recs = NULL;
		piece.keys = NULL;
		piece.n = 0;
		piece.slow = new ParallelLoad::SlowPiece();
		load.slowPieces.push_back(piece.slow);
		list.push_back(piece);
	}
	ParallelLoad::SlowPiece &slow = *list.back().slow;
	BinEntry entry;
	entry.start = start;
	entry.end = end;
	entry.line = line;
	slow.recs.push_back(entry);
	slow.keys.push_back((uint32_t)binNum | ((uint32_t)level << STAGE_LEVEL_SHIFT));
	ChromBins &bins = *_loadBins;
	if (bins.levelHi[level] < bins.levelLo[level]) {
		bins.levelLo[level] = binNum;
		bins.levelHi[level] = binNum;
	} else if (binNum < bins.levelLo[level]) {
		bins.levelLo[level] = binNum;
	} else if (binNum > bins.levelHi[level]) {
		bins.levelHi[level] = binNum;
	}
	if (entry.start == entry.end) {
		bins.coordsStable = false;
	}
}

// Indexes every chromosome loadLinesParallel loaded, largest first, on the
// pool's threads; each touches only its own ChromBins.
void BinTree::finalizeParallel()
{
	struct Job {
		ChromBins *bins;
		vector<ChromBins::Piece> pieces;
		size_t total;
	};
	ParallelLoad &load = *_parLoad;
	vector<Job> jobs(load.pieces.size());
	size_t j = 0;
	for (map<ChromBins *, vector<ParallelLoad::PieceRef> >::iterator it = load.pieces.begin();
	     it != load.pieces.end(); ++it, ++j) {
		jobs[j].bins = it->first;
		jobs[j].total = 0;
		const vector<ParallelLoad::PieceRef> &refs = it->second;
		jobs[j].pieces.resize(refs.size());
		for (size_t k = 0; k < refs.size(); k++) {
			ChromBins::Piece &piece = jobs[j].pieces[k];
			if (refs[k].slow != NULL) {
				piece.recs = refs[k].slow->recs.data();
				piece.keys = refs[k].slow->keys.data();
				piece.n = refs[k].slow->recs.size();
			} else {
				piece.recs = refs[k].recs;
				piece.keys = refs[k].keys;
				piece.n = refs[k].n;
			}
			jobs[j].total += piece.n;
		}
	}
	vector<size_t> order(jobs.size());
	for (size_t k = 0; k < order.size(); k++) {
		order[k] = k;
	}
	std::stable_sort(order.begin(), order.end(),
	                 [&jobs](size_t a, size_t b) { return jobs[a].total > jobs[b].total; });
	std::atomic<size_t> next(0);
	_pool->run([&](int) {
		for (;;) {
			size_t k = next.fetch_add(1);
			if (k >= order.size()) {
				return;
			}
			Job &job = jobs[order[k]];
			job.bins->finalizePieces(job.pieces.data(), job.pieces.size());
		}
	});
	delete _parLoad;
	_parLoad = NULL;
}

// finalize(), with the staged records in pieces.
void BinTree::ChromBins::finalizePieces(const Piece *pieces, size_t numPieces) {
	size_t numNew = 0;
	for (size_t p = 0; p < numPieces; p++) {
		numNew += pieces[p].n;
	}
	if (numNew == 0) {
		return;
	}
	vector<BinSegment> ranges;
	for (int i = 0; i < NUM_BIN_LEVELS; i++) {
		if (levelHi[i] >= levelLo[i]) {
			BinSegment r;
			r.lo = levelLo[i];
			r.hi = levelHi[i];
			r.base = 0;
			ranges.push_back(r);
		}
	}
	for (size_t a = 1; a < ranges.size(); a++) {
		BinSegment r = ranges[a];
		size_t b = a;
		while (b > 0 && ranges[b - 1].lo > r.lo) {
			ranges[b] = ranges[b - 1];
			b--;
		}
		ranges[b] = r;
	}
	segments.clear();
	for (size_t a = 0; a < ranges.size(); a++) {
		if (!segments.empty() && ranges[a].lo <= segments.back().hi + 1) {
			if (ranges[a].hi > segments.back().hi) {
				segments.back().hi = ranges[a].hi;
			}
		} else {
			segments.push_back(ranges[a]);
		}
	}
	size_t numBins = 0;
	for (size_t a = 0; a < segments.size(); a++) {
		segments[a].base = numBins;
		numBins += (size_t)(segments[a].hi - segments[a].lo + 1);
	}
	int64_t levelSlot[NUM_BIN_LEVELS];
	for (int i = 0; i < NUM_BIN_LEVELS; i++) {
		levelSlot[i] = 0;
		if (levelHi[i] < levelLo[i]) {
			continue;
		}
		for (size_t a = 0; a < segments.size(); a++) {
			if (segments[a].lo <= levelLo[i] && levelHi[i] <= segments[a].hi) {
				levelSlot[i] = (int64_t)segments[a].base - segments[a].lo;
				break;
			}
		}
	}

	const uint32_t binMask = (1u << STAGE_LEVEL_SHIFT) - 1;
	binStarts.assign(numBins + 1, 0);
	for (size_t r = 0; r < numRecs; r++) {
		int level = -1;
		binNumType bin = getBinLevel((binNumType)recs[r].start, (binNumType)recs[r].end, level);
		binStarts[(size_t)(levelSlot[level] + bin) + 1]++;
	}
	for (size_t p = 0; p < numPieces; p++) {
		const uint32_t *keys = pieces[p].keys;
		for (size_t r = 0; r < pieces[p].n; r++) {
			uint32_t key = keys[r];
			binStarts[(size_t)(levelSlot[key >> STAGE_LEVEL_SHIFT] + (int64_t)(key & binMask)) + 1]++;
		}
	}
	for (size_t k = 0; k < numBins; k++) {
		binStarts[k + 1] += binStarts[k];
	}
	vector<size_t> cursor(binStarts.begin(), binStarts.end() - 1);
	size_t total = numRecs + numNew;
	BinEntry *sorted = new BinEntry[total];
	for (size_t r = 0; r < numRecs; r++) {
		int level = -1;
		binNumType bin = getBinLevel((binNumType)recs[r].start, (binNumType)recs[r].end, level);
		sorted[cursor[(size_t)(levelSlot[level] + bin)]++] = recs[r];
	}
	for (size_t p = 0; p < numPieces; p++) {
		const uint32_t *keys = pieces[p].keys;
		const BinEntry *staged = pieces[p].recs;
		for (size_t r = 0; r < pieces[p].n; r++) {
			uint32_t key = keys[r];
			sorted[cursor[(size_t)(levelSlot[key >> STAGE_LEVEL_SHIFT] + (int64_t)(key & binMask))]++] = staged[r];
		}
	}
	delete [] recs;
	recs = sorted;
	numRecs = total;
}

// ---------------------------------------------------------------------------
// -threads: answering queries on several threads (see BinTree.h)

// getHits' Records stay unbuilt and unmoved while queries run on the
// threads: B loaded as lines is searched by coordinates alone, since its
// file has no genome file, so neither have the queries (see getHits); and
// with every chromosome's coordinates stable, printing leaves every Record's
// coordinates as they were.
bool BinTree::queriesParallel() const
{
	if (_lazy && _context->hasGenomeFile()) {
		return false;
	}
	for (mainMapType::const_iterator it = _mainMap.begin(); it != _mainMap.end(); ++it) {
		if (!it->second.coordsStable) {
			return false;
		}
	}
	return true;
}

// prepareQuery chooses BY_RECORD only for B loaded as Records with a strand
// or fraction test, or for queries with a chromosome id, which only a genome
// file gives.
bool BinTree::findReadsRecords() const
{
	return !_lazy && (_context->getSameStrand() || _context->getDiffStrand()
	                  || _context->getOverlapFractionA() != 0.0 || _context->getOverlapFractionB() != 0.0
	                  || _context->hasGenomeFile());
}

void BinTree::prepareQuery(const Record *record, QueryCursor &cursor, Query &query) const
{
	query.bins = NULL;
	query.record = record;
	query.start = record->getStartPos();
	query.end = record->getEndPos();
	query.mode = NO_HITS;
	if (record->isUnmapped()) {
		return;
	}
	if (!cursor.ready) {
		// lazily loaded B is always coordinates-only here (queriesParallel)
		cursor.allCoords = _coordHits || _lazy;
		cursor.optsCoords = !_context->getSameStrand() && !_context->getDiffStrand()
		                    && _context->getOverlapFractionA() == 0.0 && _context->getOverlapFractionB() == 0.0;
		cursor.ready = true;
	}
	const string &chr = record->getChrName();
	if (!cursor.have || chr != cursor.chr) {
		mainMapType::const_iterator mainIter = _mainMap.find(chr);
		cursor.bins = (mainIter == _mainMap.end()) ? NULL : &mainIter->second;
		cursor.chr = chr;
		cursor.have = true;
	}
	if (cursor.bins == NULL) {
		return;
	}
	query.bins = cursor.bins;
	// getHits' coordsOnly
	bool coordsOnly = cursor.allCoords
	                  || (cursor.bins->coordsStable && cursor.optsCoords && record->getChromId() == -1);
	query.mode = coordsOnly ? COORDS_ONLY : BY_RECORD;
}

// getHits' candidate runs for a query with hits, kept per thread
BinTree::QueryState &BinTree::queryRanges(const Query &query, int thread)
{
	QueryState &state = *_queryStates[thread];
	const ChromBins &bins = *query.bins;
	binNumType startBin = ((binNumType)query.start >> _binFirstShift);
	binNumType endBin = (((binNumType)query.end - 1) >> _binFirstShift);
	if (state.rangeBins != &bins || startBin != state.rangeStartBin || endBin != state.rangeEndBin) {
		const BinSegment *segBegin = bins.segments.data();
		const BinSegment *segEnd = segBegin + bins.segments.size();
		const size_t *binStarts = bins.binStarts.data();
		const BinEntry *recs = bins.recs;
		binNumType levelStart = startBin;
		binNumType levelEnd = endBin;
		int numRanges = 0;
		for (binNumType i = 0; i < NUM_BIN_LEVELS; i++) {
			binNumType offset = _binOffsetsExtended[i];
			binNumType first = levelStart + offset;
			binNumType last = levelEnd + offset;
			for (const BinSegment *seg = segBegin; seg != segEnd; ++seg) {
				if (seg->hi < first) {
					continue;
				}
				if (seg->lo > last) {
					break;
				}
				binNumType lo = first > seg->lo ? first : seg->lo;
				binNumType hi = last < seg->hi ? last : seg->hi;
				const size_t *slot = binStarts + seg->base + (size_t)(lo - seg->lo);
				if (slot[0] < slot[hi - lo + 1]) {
					state.ranges[numRanges].begin = recs + slot[0];
					state.ranges[numRanges].end = recs + slot[hi - lo + 1];
					numRanges++;
				}
			}
			levelStart >>= _binNextShift;
			levelEnd >>= _binNextShift;
		}
		state.numRanges = numRanges;
		state.rangeBins = &bins;
		state.rangeStartBin = startBin;
		state.rangeEndBin = endBin;
	}
	return state;
}

void BinTree::findHits(const Query &query, int thread, vector<BinEntry> &hits)
{
	if (query.mode == NO_HITS) {
		return;
	}
	const QueryState &state = queryRanges(query, thread);
	binNumType startPos = query.start;
	binNumType endPos = query.end;
	const RecRange *rangeEnd = state.ranges + state.numRanges;
	if (query.mode == COORDS_ONLY) {
		bool localZeroLen = (endPos - startPos == 0);
		for (const RecRange *range = state.ranges; range != rangeEnd; ++range) {
			for (const BinEntry *dbEntry = range->begin, *runEnd = range->end; dbEntry != runEnd; ++dbEntry) {
				CHRPOS maxStart = max((CHRPOS)startPos, dbEntry->start);
				CHRPOS minEnd = min((CHRPOS)endPos, dbEntry->end);
				if (minEnd < maxStart) {
					continue;
				}
				if (minEnd == maxStart && !localZeroLen && dbEntry->end - dbEntry->start != 0) {
					continue;
				}
				hits.push_back(*dbEntry);
			}
		}
		return;
	}
	// B loaded as Records
	bool sameStrand = _context->getSameStrand();
	bool diffStrand = _context->getDiffStrand();
	float overlapFractionA = _context->getOverlapFractionA();
	float overlapFractionB = _context->getOverlapFractionB();
	for (const RecRange *range = state.ranges; range != rangeEnd; ++range) {
		for (const BinEntry *dbEntry = range->begin, *runEnd = range->end; dbEntry != runEnd; ++dbEntry) {
			if (query.record->intersects(dbEntry->rec,
			                             sameStrand,
			                             diffStrand,
			                             overlapFractionA,
			                             overlapFractionB,
			                             _context->getReciprocalFraction(),
			                             _context->getEitherFraction(),
			                             _context->getObeySplits())) {
				hits.push_back(*dbEntry);
			}
		}
	}
}

// getHits' handling of the hits it finds, on the main thread.
void BinTree::serveHits(const BinEntry *hits, size_t numHits, RecordKeyVector &hitSet)
{
	ServeState &serve = *_serve;
	serve.coordUsed = 0;
	serve.lineUsed = 0;
	for (size_t k = 0; k < numHits; k++) {
		const BinEntry *dbEntry = hits + k;
		if (_coordHits) {
			CoordRecord *hit;
			if (serve.coordUsed < serve.coordPool.size()) {
				hit = serve.coordPool[serve.coordUsed];
			} else {
				hit = new CoordRecord();
				serve.coordPool.push_back(hit);
			}
			serve.coordUsed++;
			hit->setCoords(dbEntry->start, dbEntry->end);
			hitSet.push_back(hit);
			continue;
		}
		if (_lineHits) {
			LineHit *hit;
			if (serve.lineUsed < serve.linePool.size()) {
				hit = serve.linePool[serve.lineUsed];
			} else {
				hit = new LineHit(this);
				serve.linePool.push_back(hit);
			}
			serve.lineUsed++;
			hit->set(*dbEntry, _lineFileIdx[dbEntry->line->fileNum]);
			hitSet.push_back(hit);
			continue;
		}
		hitSet.push_back(_lazy ? lineRecord(dbEntry->line) : dbEntry->rec);
	}
	if (_context->getSortOutput()) {
		hitSet.sortVector();
	}
}

bool BinTree::hitsCounted() const
{
	return _coordHits && (_context->getAnyHit() || _context->getNoHit() || _context->getWriteCount());
}

// findHits' coordinates-only test, counting (_coordHits is always so).
size_t BinTree::countHits(const Query &query, int thread)
{
	if (query.mode == NO_HITS) {
		return 0;
	}
	const QueryState &state = queryRanges(query, thread);
	CHRPOS startPos = query.start;
	CHRPOS endPos = query.end;
	bool localZeroLen = (endPos - startPos == 0);
	size_t numHits = 0;
	const RecRange *rangeEnd = state.ranges + state.numRanges;
	for (const RecRange *range = state.ranges; range != rangeEnd; ++range) {
		for (const BinEntry *dbEntry = range->begin, *runEnd = range->end; dbEntry != runEnd; ++dbEntry) {
			CHRPOS maxStart = max(startPos, dbEntry->start);
			CHRPOS minEnd = min(endPos, dbEntry->end);
			if (minEnd < maxStart) {
				continue;
			}
			if (minEnd == maxStart && !localZeroLen && dbEntry->end - dbEntry->start != 0) {
				continue;
			}
			numHits++;
		}
	}
	return numHits;
}

void BinTree::serveCount(size_t numHits, RecordKeyVector &hitSet)
{
	ServeState &serve = *_serve;
	if (numHits == 0) {
		return;
	}
	if (serve.coordPool.empty()) {
		serve.coordPool.push_back(new CoordRecord());
	}
	CoordRecord *hit = serve.coordPool[0];
	for (size_t k = 0; k < numHits; k++) {
		hitSet.push_back(hit);
	}
}
