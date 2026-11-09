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
#include "Record.h"
#include "Bed3Interval.h"

using namespace std;

class FileRecordMgr;
class Record;

// A database hit rebuilt from a compact entry. With stored text (a 4-byte length, then
// the bytes of the record as it prints) print() appends that text.
class CompactHit : public Bed3Interval {
public:
	CompactHit() : _text(NULL) {}
	void setText(const char *text) { _text = text; }
	virtual void print(string &outBuf) const {
		if (_text == NULL) {
			Bed3Interval::print(outBuf);
		} else {
			uint32_t len;
			memcpy(&len, _text, sizeof(len));
			outBuf.append(_text + sizeof(len), len);
		}
	}
private:
	const char *_text;
};

typedef int64_t binNumType;
static const binNumType _binOffsetsExtended[] = {262144+32678+4096+512+64+8+1, 32678+4096+512+64+8+1, 4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};

class BinTree {
public:
	BinTree(ContextIntersect *context);

	~BinTree();
	void loadDB();
	void getHits(Record *record, RecordKeyVector &hitSet);
	// With -split and no overlap fractions, getHits drops a compact hit overlapping none of these blocks
	// (findBlockedOverlaps would drop it); NULL turns the filter off.
	void setKeyBlocks(const vector<pair<CHRPOS, CHRPOS> > *blocks) { _keyBlocks = blocks; }

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


	// A database record beside its coordinates (as stored when it was added), so a
	// candidate that cannot overlap the query is rejected without touching the Record.
	// When the database is held compactly (see loadDB) the entry carries no Record: `meta`
	// holds the file index and the zero-length flag, and getHits rebuilds a transient hit.
	struct binEntry {
		CHRPOS start;
		CHRPOS end;
		union {
			Record *record;
			int64_t meta;
		};
	};
	// For each bin number, its entries, in two parts addressed through a two-level
	// paged table (one chunk per CHUNK_BINS consecutive bin numbers, allocated on first
	// use). Part `_b` holds the bins from FIRST_BIN up (the finest level, which start-ordered
	// input fills in non-decreasing bin order): addRecordToTree appends those entries
	// straight into the final array and fills the per-bin offsets as the bin number advances.
	// Part `_a` holds the other bins: stage() collects (bin, entry) pairs and finalize() lays
	// them out in one contiguous array (stable by bin, so insertion order within a bin is
	// kept). A bin number below the last direct one converts the direct entries back to
	// staged ones and the chromosome is then laid out as one part (`_split` false).
	class allBinsType {
	public:
		ChrRef chr; // the chromosome's name, for hits rebuilt from compact entries
		allBinsType() : _split(true), _dBin(-1), _dChunk(NULL), _minBin(FIRST_BIN) {
			_a.table.assign((FIRST_BIN + CHUNK_BINS - 1) / CHUNK_BINS, (Chunk *)NULL);
			_b.table.assign((NUM_BINS - FIRST_BIN + CHUNK_BINS - 1) / CHUNK_BINS, (Chunk *)NULL);
		}
		~allBinsType() {
			for (size_t i = 0; i < _a.table.size(); i++) delete _a.table[i];
			for (size_t i = 0; i < _b.table.size(); i++) delete _b.table[i];
		}
		void stage(binNumType binNum, const binEntry &entry) {
			if (_split) {
				if (binNum < FIRST_BIN) {
					if (binNum < _minBin) _minBin = binNum;
					_stage.push_back(Staged(binNum, entry));
					return;
				}
				binNumType b = binNum - FIRST_BIN;
				if (b >= _dBin) {
					appendDirect(b, entry);
					return;
				}
				toStaged();
			}
			if (binNum < _minBin) _minBin = binNum;
			_stage.push_back(Staged(binNum, entry));
		}
		// No stored entry has a bin number below this (all of them are at or above FIRST_BIN
		// when it is FIRST_BIN), so a level whose highest reachable bin number is below it
		// cannot hold a hit.
		binNumType minBin() const { return _minBin; }
		void finalize() {
			if (_split && _dBin >= 0) closeChunk(_b.entries.size());
			size_t n = _stage.size();
			for (size_t i = 0; i < n; i++) {
				binNumType b = _stage[i].bin;
				Chunk *&chunk = _a.table[b >> CHUNK_SHIFT];
				if (chunk == NULL) chunk = new Chunk();
				chunk->off[b & (CHUNK_BINS - 1)]++;
			}
			size_t base = 0;
			for (size_t c = 0; c < _a.table.size(); c++) {
				Chunk *chunk = _a.table[c];
				if (chunk == NULL) continue;
				chunk->base = base;
				uint32_t sum = 0;
				for (int b = 0; b < CHUNK_BINS; b++) {
					uint32_t cnt = chunk->off[b];
					chunk->off[b] = sum;
					sum += cnt;
				}
				base += sum;
			}
			_a.entries.resize(n);
			// a bin is sorted (its bit is set when the chunk is made) until an entry
			// starts before the entry written just before it in the same bin
			for (size_t i = 0; i < n; i++) {
				binNumType b = _stage[i].bin;
				Chunk *chunk = _a.table[b >> CHUNK_SHIFT];
				int bb = (int)(b & (CHUNK_BINS - 1));
				uint32_t pos = chunk->off[bb]++;
				binEntry *dst = &_a.entries[chunk->base + pos];
				*dst = _stage[i].entry;
				if (pos != 0 && dst->start < dst[-1].start) chunk->sorted[bb >> 6] &= ~((uint64_t)1 << (bb & 63));
			}
			vector<Staged>().swap(_stage);
		}
		// false when the bin was never filled or is out of range
		bool find(binNumType binNum, const binEntry *&first, const binEntry *&last, bool &sorted) const {
			if ((uint64_t)binNum >= (uint64_t)NUM_BINS) return false;
			const Part *part = &_a;
			if (_split && binNum >= FIRST_BIN) {
				part = &_b;
				binNum -= FIRST_BIN;
			}
			const Chunk *chunk = part->table[binNum >> CHUNK_SHIFT];
			if (chunk == NULL) return false;
			int b = (int)(binNum & (CHUNK_BINS - 1));
			// after finalize(), off[b] is the end of bin b
			uint32_t lo = b == 0 ? 0 : chunk->off[b - 1];
			uint32_t hi = chunk->off[b];
			if (lo == hi) return false;
			first = &part->entries[chunk->base + lo];
			last = first + (hi - lo);
			sorted = (chunk->sorted[b >> 6] >> (b & 63)) & 1;
			return true;
		}
	private:
		allBinsType(const allBinsType &);
		allBinsType &operator=(const allBinsType &);
		static const int CHUNK_SHIFT = 9;
		static const binNumType CHUNK_BINS = (binNumType)1 << CHUNK_SHIFT;
		static const binNumType FIRST_BIN = 262144+32678+4096+512+64+8+1; // _binOffsetsExtended[0]
		struct Chunk {
			Chunk() : base(0) { memset(off, 0, sizeof(off)); memset(sorted, 0xff, sizeof(sorted)); }
			size_t base;
			uint64_t sorted[CHUNK_BINS / 64];
			uint32_t off[CHUNK_BINS];
		};
		struct Staged {
			Staged(binNumType b, const binEntry &e) : bin(b), entry(e) {}
			binNumType bin;
			binEntry entry;
		};
		struct Part {
			vector<Chunk *> table;
			vector<binEntry> entries;
		};
		// Appends an entry of bin number FIRST_BIN + b (b not below the last direct bin).
		void appendDirect(binNumType b, const binEntry &entry) {
			size_t n = _b.entries.size();
			if (b != _dBin) {
				if (_dBin >= 0 && (_dBin >> CHUNK_SHIFT) == (b >> CHUNK_SHIFT)) {
					// the bins between the last one and this one end where the last one ended
					uint32_t end = (uint32_t)(n - _dChunk->base);
					for (binNumType i = _dBin & (CHUNK_BINS - 1); i < (b & (CHUNK_BINS - 1)); i++) _dChunk->off[i] = end;
				} else {
					if (_dBin >= 0) closeChunk(n);
					_dChunk = new Chunk();
					_dChunk->base = n;
					_b.table[b >> CHUNK_SHIFT] = _dChunk;
				}
				_dBin = b;
			} else if (entry.start < _b.entries[n - 1].start) {
				int bb = (int)(b & (CHUNK_BINS - 1));
				_dChunk->sorted[bb >> 6] &= ~((uint64_t)1 << (bb & 63));
			}
			_b.entries.push_back(entry);
		}
		// The open bin and every later bin of its chunk end at n.
		void closeChunk(size_t n) {
			uint32_t end = (uint32_t)(n - _dChunk->base);
			for (binNumType i = _dBin & (CHUNK_BINS - 1); i < CHUNK_BINS; i++) _dChunk->off[i] = end;
		}
		// Gives up the direct part: its entries join the staged ones in arrival order
		// (their bins are above every staged bin, so per-bin order is unchanged).
		void toStaged() {
			_stage.reserve(_stage.size() + _b.entries.size());
			for (size_t i = 0; i < _b.entries.size(); i++) {
				const binEntry &e = _b.entries[i];
				_stage.push_back(Staged(BinTree::getBin(e.start, e.end), e));
			}
			vector<binEntry>().swap(_b.entries);
			for (size_t i = 0; i < _b.table.size(); i++) delete _b.table[i];
			vector<Chunk *>().swap(_b.table);
			_a.table.resize((NUM_BINS + CHUNK_BINS - 1) / CHUNK_BINS, (Chunk *)NULL);
			_dChunk = NULL;
			_split = false;
		}
		bool _split;
		binNumType _dBin; // bin (minus FIRST_BIN) of the last direct entry, -1 for none
		Chunk *_dChunk; // the chunk of that bin
		binNumType _minBin; // lowest bin number of any entry staged or placed
		Part _a, _b;
		vector<Staged> _stage;
	};
	typedef map<string, allBinsType *> mainMapType; //for each chrom, a map of bin num to RecordLists.
	mainMapType _mainMap;

	// last chrom looked up (queries) / inserted (load), to skip the map find.
	ChrRef _lastQueryChr;
	allBinsType *_lastQueryBins;
	bool _haveQueryChr;
	// per level: first entry of the last visited start bin and how many leading entries
	// had end < the previous query's start (valid while query starts do not decrease)
	const binEntry *_curFirst[8];
	size_t _curCnt[8];
	binNumType _curStart;
	// per level: the last start-bin lookup (bin number, find() result); -1 = none
	binNumType _memoBin[8];
	bool _memoFound[8];
	bool _memoSorted[8];
	const binEntry *_memoFirst[8];
	const binEntry *_memoLast[8];
	string _lastAddChr;
	allBinsType *_lastAddBins;

	map<binNumType, int> _binsHit;

	// Compact mode: the database is staged as (start, end, file index, zero-length flag)
	// and each query's hits are transient Records from a pool reset by the next getHits.
	bool _compact;
	vector<CompactHit *> _pool;
	vector<CompactHit *> _poolBlocks;
	size_t _poolUsed;
	CompactHit *nextPoolRecord();

	// Text mode (compact, and the database record is printed): each line's printed text
	// is kept in an arena of fixed-size blocks; an entry's `meta` holds
	// (offset << 16) | (fileIdx << 1) | zero-length flag.
	bool _text;
	static const size_t TEXT_BLOCK = (size_t)1 << 22;
	vector<char *> _textBlocks;
	size_t _textUsed; // bytes used in the last block
	string _textScratch;
	int64_t appendText(Record *record);
	int64_t appendRaw(const char *text, size_t len);
	int64_t _lastTextOff; // offset of the last appended text
	const vector<pair<CHRPOS, CHRPOS> > *_keyBlocks;
	bool _stopAtFirst; // only whether a query has a hit is read (-u, -v)
	bool textEquals(int64_t off, const char *text, size_t len) const;

	bool addRecordToTree(Record *);
	bool stageCoords(binNumType startPos, binNumType endPos, int64_t meta);
	static binNumType getBin(binNumType start, binNumType end);
	binNumType getBin(const Record *record) const;


};


#endif /* BINTREE_H_ */
