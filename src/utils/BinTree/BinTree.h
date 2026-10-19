/*
 * BinTree.h
 *
 *  Created on: Jan 5, 2013
 *      Author: nek3d
 */

#ifndef BINTREE_H_
#define BINTREE_H_

#include "BedtoolsTypes.h"

#include <stdint.h>
#include <string>
#include <set>
#include <map>
#include <vector>

#include "string.h"
#include "RecordKeyList.h"
#include "ContextIntersect.h"

using namespace std;

class FileRecordMgr;
class Record;
class WorkerPool;

typedef int64_t binNumType;
static const binNumType _binOffsetsExtended[] = {262144+32678+4096+512+64+8+1, 32678+4096+512+64+8+1, 4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};

class BinTree {
	struct ChromBins;
	struct LazyLine;
public:
	BinTree(ContextIntersect *context);

	~BinTree();
	void loadDB();
	void getHits(Record *record, RecordKeyVector &hitSet);

	// A record and its coordinates as they were when it was loaded, so that
	// getHits can test overlap without reading the record itself.
	struct BinEntry {
		CHRPOS start;
		CHRPOS end;
		union {
			Record *rec;    // loaded eagerly (!_lazy)
			LazyLine *line; // loaded lazily (_lazy); NULL with _coordHits
		};
	};

	// -threads N (IntersectFile). Given the pool before loadDB, loadDB parses
	// B text files on its threads (loadLinesParallel) and indexes the
	// chromosomes on them.
	void setPool(WorkerPool *pool);

	// Queries answered on the pool's threads once loadDB is done, which
	// queriesParallel allows. On the main thread, in input order,
	// prepareQuery looks up a query's chromosome with that thread's cursor
	// and records what any thread needs to search it, so that no thread but
	// the main one reads the query's Record on the common path. findHits on
	// pool thread t then appends the entries of the hits getHits would find,
	// in the same order, reading only the query's description, what loadDB
	// built and, for B loaded as Records with a strand or fraction test,
	// the Records, none of which changes while it runs. serveHits, on the
	// main thread, gives hitSet what getHits would have for the entries.
	// The main thread's: the last query's chromosome and its bins, and the
	// tree's settings prepareQuery reads, copied on its first call, so that
	// it reads nothing per query that a thread printing writes.
	struct QueryCursor {
		QueryCursor() : bins(NULL), have(false), ready(false), allCoords(false), optsCoords(false) {}
		string chr;
		const ChromBins *bins;
		bool have;
		bool ready;
		bool allCoords;  // every query coordinates-only (_coordHits or _lazy)
		bool optsCoords; // no strand or fraction test
	};
	enum QueryMode { NO_HITS, COORDS_ONLY, BY_RECORD };
	struct Query {
		const ChromBins *bins;
		const Record *record;
		CHRPOS start;
		CHRPOS end;
		QueryMode mode;
	};
	bool queriesParallel() const;
	// Whether findHits can read B's Records (BY_RECORD), which printing
	// changes for a moment (the zero-length adjustment, undone and redone),
	// so that no query may be searched while another prints.
	bool findReadsRecords() const;
	void prepareQuery(const Record *record, QueryCursor &cursor, Query &query) const;
	void findHits(const Query &query, int thread, vector<BinEntry> &hits);
	void serveHits(const BinEntry *hits, size_t numHits, RecordKeyVector &hitSet);
	// Where intersect reads only how many hits a query has (-u, -v, -c with
	// coordinate stand-ins), the threads count them (countHits) and
	// serveCount gives hitSet that many stand-ins, all the same one.
	bool hitsCounted() const;
	size_t countHits(const Query &query, int thread);
	void serveCount(size_t numHits, RecordKeyVector &hitSet);

private:

	ContextIntersect *_context;

    //
    // BIN HANDLING
    //
        typedef int64_t binNumType;

#define PRId_BINNUMTYPE PRId64

	// bins range in size from 16kb to 32Gb
	static const binNumType NUM_BIN_LEVELS = 8;

	static const binNumType _binFirstShift = 14;       /* How much to shift to get to finest bin. */
	static const binNumType _binNextShift  = 3;        /* How much to shift to get to next larger bin. */

	static const binNumType NUM_BINS = (1 << (_binNextShift * NUM_BIN_LEVELS)) / ((1 << _binNextShift) - 1);


	// A run of consecutive global bin numbers [lo, hi] with one slot each in
	// ChromBins::binStarts, starting at base.
	struct BinSegment {
		binNumType lo;
		binNumType hi;
		size_t base;
	};

	// A B record loaded without its Record (see loadDB): the line's text
	// follows this header in the arena, and rec is built from it when the
	// record is first a hit, or at once if a warning has to print it.
	// zeroLength and verbatim are set only for line hits (_lineHits). With
	// coordinate stand-ins (_coordHits) no line is kept.
	struct LazyLine {
		Record *rec;
		uint32_t len;
		uint32_t fileNum : 30;   // its database file
		uint32_t zeroLength : 1; // loaded with start == end
		uint32_t verbatim : 1;   // its Record, adjustment undone, prints as the text
		const char *text() const { return reinterpret_cast<const char *>(this + 1); }
	};

	// Records staged for indexing and their keys (bin | (level <<
	// STAGE_LEVEL_SHIFT)), in insertion order, in blocks that are filled once
	// and never moved. A block of size class c holds STAGE_BLOCK_MIN << c
	// entries; a chromosome's blocks double in size up to the largest class.
	// The blocks a chromosome has indexed go to the tree's pool, from which
	// later chromosomes take them, so that staging reuses the same memory.
	static const size_t STAGE_BLOCK_MIN = 8;
	static const int STAGE_SIZE_CLASSES = 14; // up to 65536 entries
	struct StageBlock {
		BinEntry *recs;
		uint32_t *keys;
		int sizeClass;
	};
	struct StagePool {
		vector<StageBlock> blocks[STAGE_SIZE_CLASSES]; // free blocks by size class
		StagePool() {}
		~StagePool() { release(); }
		StageBlock take(int sizeClass);
		void give(const StageBlock &block) { blocks[block.sizeClass].push_back(block); }
		void release();
		StagePool(const StagePool &) = delete;
		StagePool &operator=(const StagePool &) = delete;
	};
	StagePool _stagePool;

	// One chromosome's records, indexed by global bin number. While loading,
	// records are staged in insertion order with their bin and level; finalize()
	// then counting-sorts them (stably) into recs, so a bin's records are
	// contiguous and keep insertion order. The bins of segment s are
	// binStarts[s.base .. s.base + s.hi - s.lo]; bin j's records are
	// recs[binStarts[k]] to recs[binStarts[k + 1] - 1], k = s.base + j - s.lo.
	// Segments are the per-level ranges of bins used, merged where they meet,
	// so lookups go by global bin number however the levels' ranges lie.
	// finalize() can run more than once: records staged since the last one
	// follow those already in recs, which were all inserted before them.
	struct ChromBins {
		ChromBins();
		~ChromBins();
		void stage(const BinEntry &entry, binNumType binNum, int level, StagePool &pool);
		void finalize(StagePool &pool);

		// the blocks staged into, in order: all full but the last, which has
		// stageLeft slots free
		vector<StageBlock> stageBlocks;
		BinEntry *stageRecPos;
		uint32_t *stageKeyPos;
		size_t stageLeft;
		size_t numStaged;
		void addStageBlock(StagePool &pool);
		binNumType levelLo[NUM_BIN_LEVELS];
		binNumType levelHi[NUM_BIN_LEVELS];

		vector<BinSegment> segments; // sorted by lo, disjoint
		vector<size_t> binStarts;
		// numRecs entries (NULL before the first finalize), allocated
		// uninitialized, since the counting sort writes every one
		BinEntry *recs;
		size_t numRecs;

		// it owns recs and its blocks; std::map builds it in place
		ChromBins(const ChromBins &) = delete;
		ChromBins &operator=(const ChromBins &) = delete;

		// True when no record was loaded with start == end. Output undoes and
		// redoes the zero-length adjustment on every hit it prints, which
		// restores a record's coordinates unless it has start == end without
		// having been adjusted; then adjustZeroLength moves them for good. So
		// while this holds, the coordinates in recs are the records' own.
		bool coordsStable;

		// finalize() for loadLinesParallel, whose records are staged in the
		// given pieces, in insertion order, rather than in stageBlocks
		struct Piece;
		void finalizePieces(const Piece *pieces, size_t numPieces);
	};
	// Bins are below NUM_BINS < 2^22, so the level fits above them.
	static const int STAGE_LEVEL_SHIFT = 22;

	typedef map<string, ChromBins> mainMapType; //for each chrom, its bin index.
	mainMapType _mainMap;

	// The chromosome of the last record loaded and of the last query, with its
	// bin index in _mainMap (NULL for a query whose chromosome is not there).
	// Inputs usually come grouped by chromosome, so consecutive records skip
	// the map search; pointers to map values stay valid across inserts.
	string _loadChr;
	ChromBins *_loadBins;
	string _queryChr;
	const ChromBins *_queryBins;
	bool _haveQuery;

	// The last query's candidates: for each level and segment, in the order
	// getHits visits them, the non-empty run of recs holding the bins the
	// query covers there. Those bins follow from the chromosome and the
	// query's first and last level-0 bins alone, so a query with the same
	// three scans the same runs. A chromosome has at most NUM_BIN_LEVELS
	// segments, so a query has at most MAX_RANGES runs. _rangeBins is NULL
	// when none are kept.
	struct RecRange {
		const BinEntry *begin;
		const BinEntry *end;
	};
	static const int MAX_RANGES = NUM_BIN_LEVELS * NUM_BIN_LEVELS;
	const ChromBins *_rangeBins;
	binNumType _rangeStartBin;
	binNumType _rangeEndBin;
	int _numRanges;
	RecRange _ranges[MAX_RANGES];

	map<binNumType, int> _binsHit;

	// Whether B was loaded as lines (every entry's line set) or as Records
	// (every entry's rec set), and the arena holding the lines.
	bool _lazy;

	// Whether a hit of B loaded as lines can be a coordinates-only stand-in
	// in place of the line's Record (see loadDB), and the stand-ins: the
	// first _coordUsed of _coordPool are the current query's hits.
	class CoordRecord;
	bool _coordHits;
	vector<CoordRecord *> _coordPool;
	size_t _coordUsed;

	// Whether a hit of B loaded as lines can be a line stand-in in place of
	// its Record where intersect prints or tags hits (see loadDB), the
	// stand-ins (as _coordPool), and each database file's file index.
	class LineHit;
	bool _lineHits;
	vector<LineHit *> _linePool;
	size_t _lineUsed;
	vector<int> _lineFileIdx;

	vector<char *> _arenaChunks;
	char *_arenaPos;
	size_t _arenaLeft;
	static const size_t ARENA_CHUNK = 4 << 20;

	bool linesLoadable();
	void loadLines(int fileNum, FileRecordMgr *databaseFile);

	// -threads: the pool, loadLinesParallel's state while loading (NULL
	// when B is loaded on one thread), and findHits' state for each thread
	WorkerPool *_pool;
	struct ParallelLoad;
	struct LoadBlock;
	struct LoadFormat;
	ParallelLoad *_parLoad;
	struct QueryState;
	vector<QueryState *> _queryStates;
	QueryState &queryRanges(const Query &query, int thread);
	// serveHits' and serveCount's stand-ins, apart from getHits' and from
	// what the other threads read (ServeState)
	struct ServeState;
	ServeState *_serve;
	void loadLinesParallel(int fileNum, FileRecordMgr *databaseFile);
	void parseBlock(LoadBlock &block, const LoadFormat &format) const;
	void consumeBlock(LoadBlock &block, int fileNum, FileRecordMgr *databaseFile);
	void loadSlowLine(const char *text, size_t len, int fileNum, FileRecordMgr *databaseFile, bool asText);
	void setLoadChrom(const char *chr, size_t len);
	void finalizeParallel();
	LazyLine *storeLine(const char *text, size_t len, int fileNum, bool zeroLength, bool verbatim);
	Record *lineRecord(LazyLine *line) {
		return line->rec != NULL ? line->rec : buildLineRecord(line);
	}
	Record *buildLineRecord(LazyLine *line);

	bool addRecordToTree(Record *);
	bool addEntryToTree(const string &chr, const BinEntry &entry);
	binNumType getBin(binNumType start, binNumType end) const;
	binNumType getBin(const Record *record) const;
	static binNumType getBinLevel(binNumType start, binNumType end, int &level);


};


#endif /* BINTREE_H_ */
