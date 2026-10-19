/*****************************************************************************
  intersectFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "intersectFile.h"
#include "ContextIntersect.h"
#include "FileRecordMgr.h"
#include "BinTree.h"
#include "RecordOutputMgr.h"
#include "WorkerPool.h"

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <system_error>


// Unsorted queries on -threads N threads. The main thread reads the queries
// in batches exactly as nextUnsortedFind reads them one at a time, so
// records, warnings and errors come from the query file in the same order,
// and describes each for the search (BinTree::prepareQuery). The pool's
// threads, the main one included, then find a batch's hits
// (BinTree::findHits, or countHits where only counts matter), in rounds
// that stop taking queries past a budget of hits. A printer thread takes
// the rounds in input order and does for each query what the driver loop
// does (BinTree::serveHits, checkSplits, processHits) while the main thread
// reads the next batch; the printer alone prints, and the main thread alone
// allocates and frees query Records, recycling a batch's once it is printed.
// Hits are found only for queries not yet printing, from what loadDB built,
// which printing leaves as it was (BinTree::queriesParallel).
//
// An error that exits while a batch is read would lose the queries read
// before it, which at one thread were printed, the output flushed as far as
// it then was, before the bad record was read. An exit handler waits for the
// printer to print every round it was given, then prints the queries read
// since, so the output is the same bytes.
class IntersectFile::ParallelQueries {
public:
	ParallelQueries(IntersectFile *tool, WorkerPool *pool);
	~ParallelQueries();
	// Answers and prints every query on its first call, then returns false.
	bool next(RecordKeyVector &hits);
	// with key lines, the reader whose current line is that of the key the
	// printer prints
	const SingleLineDelimTextFileReader *keyLineReader() const { return _keyReader; }

private:
	static const size_t BATCH = 16384;        // queries read at once
	static const size_t CHUNK = 64;           // queries per share of the hit search
	static const size_t HIT_BUDGET = 1 << 22; // a round stops taking shares past this many hits
	static const size_t NUM_SLOTS = 4;        // batches read or printing
	static const size_t NUM_ROUNDS = 3;       // rounds found or printing
	// The memory the printer reads, the main thread writes again only
	// batches later: a line the printer still holds costs the main thread a
	// wait for the printer's copy to be dropped, which between cores that
	// share no cache doubled the time the main thread took to read. So the
	// key lines rotate through more sets than there are slots, and a
	// printed batch's Records are freed (for the free list to give out
	// again) RECORD_DELAY batches later.
	static const size_t NUM_KEY_LINE_SETS = 16;
	static const size_t RECORD_DELAY = 4;

	// The hits of up to CHUNK consecutive queries: query j's are
	// hits[begin[j]] to hits[begin[j + 1] - 1]; where only counts matter
	// (_counted), begin holds running counts and hits stays empty.
	struct Chunk {
		vector<BinTree::BinEntry> hits;
		vector<size_t> begin;
	};
	// With key lines, a batch's query lines as the query reader read them:
	// line i's bytes and a NUL from start[i], then its numFields + 1
	// delimiter positions from delims[i * (numFields + 1)]. (This, Slot and
	// Round begin with a line's worth of padding, so that the main thread
	// writing one never shares a line with the printer reading another.)
	struct KeyLines {
		char pad[64];
		vector<char> text;
		vector<size_t> start;
		vector<int> delims;
	};
	// A batch of queries.
	struct Slot {
		char pad[64];
		vector<Record *> batch;
		vector<BinTree::Query> queries;
		KeyLines *lines;
		bool busy; // read and not yet all printed
		Slot() : lines(NULL), busy(false) {}
	};
	// Queries [first, last) of a slot with their hits found.
	struct Round {
		char pad[64];
		Slot *slot;
		size_t first;
		size_t last;
		bool lastOfSlot;
		vector<Chunk> chunks;
		bool busy; // found and not yet printed
		Round() : slot(NULL), first(0), last(0), lastOfSlot(false), busy(false) {}
	};

	IntersectFile *_tool;
	BinTree *_binTree;
	FileRecordMgr *_queryFRM;
	ContextIntersect *_context;
	WorkerPool *_pool;
	bool _counted; // only hit counts matter (BinTree::hitsCounted)
	int _keyFields;
	SingleLineDelimTextFileReader *_keyReader; // with key lines, the printer's

	Slot _slots[NUM_SLOTS];
	Round _rounds[NUM_ROUNDS];
	KeyLines _keyLineSets[NUM_KEY_LINE_SETS];
	vector<Record *> _retired[RECORD_DELAY]; // printed, waiting to be freed
	BinTree::QueryCursor _cursor; // the main thread's, for prepareQuery
	Slot *_readSlot;              // the slot the main thread is reading into

	// the round being found: its queries' shares are taken from _nextChunk;
	// the counters the finding threads write are kept apart from what the
	// printer reads
	Round *_findRound;
	size_t _numChunks;
	char _pad[64];
	std::atomic<size_t> _nextChunk;
	std::atomic<size_t> _numHits;
	char _pad2[64];

	// rounds handed to the printer and printed (in order: round k is
	// _rounds[k % NUM_ROUNDS]); guarded by _mutex, as are the busy flags
	std::thread _printer;
	bool _printerRunning;
	size_t _roundsMade;
	size_t _roundsPrinted;
	bool _finished;
	std::mutex _mutex;
	std::condition_variable _toPrinter;
	std::condition_variable _toMain;

	static ParallelQueries *_reading; // reading a batch, for atExit
	static void atExit();

	void run();
	bool readBatch(Slot &slot, size_t b);
	void recycle(Slot &slot, size_t b);
	void findShare(int thread);
	void handRound(Round &round);
	void printLoop();
	void printRound(Round &round);
	void showKeyLine(const Slot &slot, size_t i);
	void finishRead();
};

IntersectFile::ParallelQueries *IntersectFile::ParallelQueries::_reading = NULL;

IntersectFile::ParallelQueries::ParallelQueries(IntersectFile *tool, WorkerPool *pool)
: _tool(tool),
  _binTree(tool->_binTree),
  _queryFRM(tool->_queryFRM),
  _context(tool->upCast(tool->_context)),
  _pool(pool),
  _counted(tool->_binTree->hitsCounted()),
  _keyFields(0),
  _keyReader(NULL),
  _readSlot(NULL),
  _findRound(NULL),
  _numChunks(0),
  _nextChunk(0),
  _numHits(0),
  _printerRunning(false),
  _roundsMade(0),
  _roundsPrinted(0),
  _finished(false)
{
	static bool handlerSet = false;
	if (!handlerSet) {
		atexit(atExit);
		handlerSet = true;
	}
	if (tool->_keyLines) {
		const SingleLineDelimTextFileReader *reader = _queryFRM->textReader();
		_keyFields = reader->getNumFields();
		_keyReader = new SingleLineDelimTextFileReader(_keyFields, reader->getDelimChar());
	}
}

IntersectFile::ParallelQueries::~ParallelQueries()
{
	delete _keyReader;
}

bool IntersectFile::ParallelQueries::next(RecordKeyVector &)
{
	if (!_finished) {
		run();
	}
	return false;
}

void IntersectFile::ParallelQueries::run()
{
	if (_tool->_outputMgr != NULL && !_binTree->findReadsRecords()) {
		try {
			_printer = std::thread(&ParallelQueries::printLoop, this);
			_printerRunning = true;
		} catch (const std::system_error &) {
			// rounds print on this thread as they are found
		}
	}
	for (size_t b = 0; ; b++) {
		Slot &slot = _slots[b % NUM_SLOTS];
		{
			std::unique_lock<std::mutex> lock(_mutex);
			while (slot.busy) {
				_toMain.wait(lock);
			}
		}
		recycle(slot, b);
		if (!readBatch(slot, b)) {
			break;
		}
		{
			std::lock_guard<std::mutex> lock(_mutex);
			slot.busy = true;
		}
		size_t ready = 0;
		while (ready < slot.batch.size()) {
			Round &round = _rounds[_roundsMade % NUM_ROUNDS];
			{
				std::unique_lock<std::mutex> lock(_mutex);
				while (round.busy) {
					_toMain.wait(lock);
				}
			}
			// the hits of the slot's queries from ready on, as far as
			// HIT_BUDGET lets the threads go: shares are taken in order, so
			// the queries done are those of the first shares
			round.slot = &slot;
			round.first = ready;
			_numChunks = (slot.batch.size() - ready + CHUNK - 1) / CHUNK;
			if (round.chunks.size() < _numChunks) {
				round.chunks.resize(_numChunks);
			}
			_findRound = &round;
			_nextChunk = 0;
			_numHits = 0;
			_pool->run([this](int thread) { findShare(thread); });
			size_t chunksDone = std::min(_nextChunk.load(), _numChunks);
			ready = std::min(ready + chunksDone * CHUNK, slot.batch.size());
			round.last = ready;
			round.lastOfSlot = (ready == slot.batch.size());
			handRound(round);
		}
	}
	{
		std::lock_guard<std::mutex> lock(_mutex);
		_finished = true;
	}
	_toPrinter.notify_all();
	if (_printerRunning) {
		_printer.join();
		_printerRunning = false;
	}
	for (size_t s = 0; s < NUM_SLOTS + RECORD_DELAY; s++) {
		recycle(_slots[s % NUM_SLOTS], s);
	}
}

// cleanupHits' frees, for the Records of the batch the printer finished
// RECORD_DELAY batches before batch b's slot was last read into
void IntersectFile::ParallelQueries::recycle(Slot &slot, size_t b)
{
	vector<Record *> &retired = _retired[b % RECORD_DELAY];
	for (size_t i = 0; i < retired.size(); i++) {
		_queryFRM->deleteRecord(retired[i]);
	}
	retired.swap(slot.batch);
	slot.batch.clear();
}

// nextUnsortedFind's reading, for up to BATCH queries
bool IntersectFile::ParallelQueries::readBatch(Slot &slot, size_t b)
{
	// batch b - NUM_KEY_LINE_SETS, which last used these, is printed, as
	// batch b - NUM_SLOTS is
	KeyLines &lines = _keyLineSets[b % NUM_KEY_LINE_SETS];
	lines.text.clear();
	lines.start.assign(1, 0);
	lines.delims.clear();
	slot.lines = &lines;
	_readSlot = &slot;
	_reading = this;
	while (slot.batch.size() < BATCH && !_queryFRM->eof()) {
		Record *queryRecord = _queryFRM->getNextRecord();
		if (queryRecord == NULL) {
			continue;
		}
		_context->testNameConventions(queryRecord);
		slot.batch.push_back(queryRecord);
		if (_keyReader != NULL) {
			const SingleLineDelimTextFileReader *reader = _queryFRM->textReader();
			const char *line = reader->currentLine();
			size_t len = reader->currentLineLength();
			lines.text.insert(lines.text.end(), line, line + len + 1);
			lines.start.push_back(lines.text.size());
			const int *delims = reader->delimPositions();
			lines.delims.insert(lines.delims.end(), delims, delims + _keyFields + 1);
		}
	}
	_reading = NULL;
	slot.queries.resize(slot.batch.size());
	for (size_t i = 0; i < slot.batch.size(); i++) {
		_binTree->prepareQuery(slot.batch[i], _cursor, slot.queries[i]);
	}
	return !slot.batch.empty();
}

void IntersectFile::ParallelQueries::findShare(int thread)
{
	Round &round = *_findRound;
	const Slot &slot = *round.slot;
	for (;;) {
		if (_numHits.load(std::memory_order_relaxed) >= HIT_BUDGET) {
			return;
		}
		size_t c = _nextChunk.fetch_add(1);
		if (c >= _numChunks) {
			return;
		}
		Chunk &chunk = round.chunks[c];
		chunk.hits.clear();
		chunk.begin.clear();
		size_t first = round.first + c * CHUNK;
		size_t last = std::min(first + CHUNK, slot.batch.size());
		if (_counted) {
			size_t count = 0;
			for (size_t q = first; q < last; q++) {
				chunk.begin.push_back(count);
				count += _binTree->countHits(slot.queries[q], thread);
			}
			chunk.begin.push_back(count);
			continue;
		}
		for (size_t q = first; q < last; q++) {
			chunk.begin.push_back(chunk.hits.size());
			_binTree->findHits(slot.queries[q], thread, chunk.hits);
		}
		chunk.begin.push_back(chunk.hits.size());
		_numHits.fetch_add(chunk.hits.size(), std::memory_order_relaxed);
	}
}

void IntersectFile::ParallelQueries::handRound(Round &round)
{
	if (!_printerRunning) {
		printRound(round);
		if (round.lastOfSlot) {
			round.slot->busy = false;
		}
		_roundsMade++;
		_roundsPrinted++;
		return;
	}
	{
		std::lock_guard<std::mutex> lock(_mutex);
		round.busy = true;
		_roundsMade++;
	}
	_toPrinter.notify_one();
}

void IntersectFile::ParallelQueries::printLoop()
{
	for (size_t k = 0; ; k++) {
		Round *round;
		{
			std::unique_lock<std::mutex> lock(_mutex);
			while (k == _roundsMade && !_finished) {
				_toPrinter.wait(lock);
			}
			if (k == _roundsMade) {
				return;
			}
			round = &_rounds[k % NUM_ROUNDS];
		}
		printRound(*round);
		{
			std::lock_guard<std::mutex> lock(_mutex);
			round->busy = false;
			if (round->lastOfSlot) {
				round->slot->busy = false;
			}
			_roundsPrinted = k + 1;
		}
		_toMain.notify_all();
	}
}

// What the driver loop does for each of the round's queries, but the free.
void IntersectFile::ParallelQueries::printRound(Round &round)
{
	const Slot &slot = *round.slot;
	RecordKeyVector hits;
	for (size_t i = round.first; i < round.last; i++) {
		hits.setKey(slot.batch[i]);
		showKeyLine(slot, i);
		const Chunk &chunk = round.chunks[(i - round.first) / CHUNK];
		size_t j = (i - round.first) % CHUNK;
		if (_counted) {
			_binTree->serveCount(chunk.begin[j + 1] - chunk.begin[j], hits);
		} else {
			_binTree->serveHits(chunk.hits.data() + chunk.begin[j], chunk.begin[j + 1] - chunk.begin[j], hits);
		}
		_tool->checkSplits(hits);
		_tool->processHits(_tool->_outputMgr, hits);
		hits.clearAll();
	}
}

void IntersectFile::ParallelQueries::showKeyLine(const Slot &slot, size_t i)
{
	if (_keyReader != NULL) {
		const KeyLines &lines = *slot.lines;
		size_t from = lines.start[i];
		_keyReader->showLine(&lines.text[from], lines.start[i + 1] - from - 1,
		                     &lines.delims[i * (_keyFields + 1)]);
	}
}

// At exit while reading a batch, on the main thread: once the printer has
// printed every round it was given, what the driver loop does for each
// query read since.
void IntersectFile::ParallelQueries::finishRead()
{
	if (_tool->_outputMgr == NULL) {
		return;
	}
	if (_printerRunning) {
		std::unique_lock<std::mutex> lock(_mutex);
		while (_roundsPrinted < _roundsMade) {
			_toMain.wait(lock);
		}
	}
	Slot &slot = *_readSlot;
	vector<BinTree::BinEntry> found;
	RecordKeyVector hits;
	for (size_t i = 0; i < slot.batch.size(); i++) {
		hits.setKey(slot.batch[i]);
		showKeyLine(slot, i);
		BinTree::Query query;
		_binTree->prepareQuery(slot.batch[i], _cursor, query);
		found.clear();
		_binTree->findHits(query, 0, found);
		_binTree->serveHits(found.data(), found.size(), hits);
		_tool->checkSplits(hits);
		_tool->processHits(_tool->_outputMgr, hits);
		hits.clearAll();
	}
}

void IntersectFile::ParallelQueries::atExit()
{
	ParallelQueries *reading = _reading;
	if (reading != NULL) {
		_reading = NULL;
		reading->finishRead();
	}
}


IntersectFile::IntersectFile(ContextIntersect *context)
: ToolBase(upCast(context)),
  _sweep(NULL),
  _binTree(NULL),
  _queryFRM(NULL),
  _keyLines(false),
  _keyLineReader(NULL),
  _outputMgr(NULL),
  _pool(NULL),
  _parallel(NULL)
{

}

IntersectFile::~IntersectFile(void) {
	delete _parallel;
	_parallel = NULL;

	delete _sweep;
	_sweep = NULL;

	delete _binTree;
	_binTree = NULL;

	delete _pool;
	_pool = NULL;
}

bool IntersectFile::init() {

	_queryFRM = upCast(_context)->getFile(upCast(_context)->getQueryFileIdx());

	 if (upCast(_context)->getSortedInput()) {
		 makeSweep();
		return _sweep->init();
	 } else {
		_binTree = new BinTree( upCast(_context));
		if (upCast(_context)->getNumThreads() > 1) {
			_pool = new WorkerPool(upCast(_context)->getNumThreads());
			_binTree->setPool(_pool);
		}
		_binTree->loadDB();
		// Unsorted, nothing reads the query file between findNext and
		// processHits, so a key's line is still the reader's current line
		// when it prints; a file linesPrintAsText allows prints a Record as
		// its fields in file order joined by tabs. Queries read in batches
		// keep their lines (ParallelQueries).
		_keyLines = _context->getProgram() == ContextBase::INTERSECT
		            && _queryFRM->linesLoadable() && _queryFRM->linesPrintAsText();
		_keyLineReader = _keyLines ? _queryFRM->textReader() : NULL;
		if (_pool != NULL && _binTree->queriesParallel()) {
			_parallel = new ParallelQueries(this, _pool);
			_keyLineReader = _parallel->keyLineReader();
		}
	 }

	 return true;
}

bool IntersectFile::findNext(RecordKeyVector &hits)
{
	 bool retVal = false;
	 if (upCast(_context)->getSortedInput()) {
		retVal = nextSortedFind(hits);
	 } 
	 else if (_parallel != NULL) {
		retVal = _parallel->next(hits);
	 }
	 else {
		retVal = nextUnsortedFind(hits);
	 }

	 if (retVal) {
		 checkSplits(hits);
	 }
	 return retVal;
}

void IntersectFile::processHits(RecordOutputMgr *outputMgr, RecordKeyVector &hits)
{
	if (_keyLines) {
		outputMgr->setKeyLine(hits.getKey(), _keyLineReader);
		outputMgr->printRecord(hits);
		outputMgr->clearKeyLine();
		return;
	}
	outputMgr->printRecord(hits);
}

void IntersectFile::cleanupHits(RecordKeyVector &hits)
{
	_queryFRM->deleteRecord(hits.getKey());
	hits.clearAll();
}

bool IntersectFile::finalizeCalculations()
{
    if (upCast(_context)->getSortedInput() && !upCast(_context)->hasGenomeFile()) 
    {
        if (_context->getNameCheckDisabled())
            _sweep->closeOut(false);
        else
            _sweep->closeOut(true);
    }
    return true;
}

bool IntersectFile::nextSortedFind(RecordKeyVector &hits)
{
    if (!_sweep->next(hits)) {
    	return false;
    }
    return true;
}

bool IntersectFile::nextUnsortedFind(RecordKeyVector &hits)
{

	while (!_queryFRM->eof()) {
		Record *queryRecord = _queryFRM->getNextRecord();
		if (queryRecord == NULL) {
			continue;
		} else {
			_context->testNameConventions(queryRecord);
			hits.setKey(queryRecord);
			 _binTree->getHits(queryRecord, hits);
			return true;
		}
	}
	return false;
}

void IntersectFile::makeSweep() {
	_sweep = new NewChromSweep(upCast(_context));
}

void IntersectFile::checkSplits(RecordKeyVector &hitSet)
{
	if (upCast(_context)->getObeySplits()) {
		RecordKeyVector keySet(hitSet.getKey());
		RecordKeyVector resultSet(hitSet.getKey());
		upCast(_context)->getSplitBlockInfo()->findBlockedOverlaps(keySet, hitSet, resultSet);
		hitSet.swap(resultSet);
	}
}
