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
	// For each bin number, its entries. Filled in two steps: addRecordToTree stages
	// (bin, entry) pairs; finalize() lays the entries of every bin out in one
	// contiguous array (stable by bin, so insertion order within a bin is kept),
	// addressed through a two-level paged table: one chunk per CHUNK_BINS
	// consecutive bin numbers, allocated on first use.
	class allBinsType {
	public:
		ChrRef chr; // the chromosome's name, for hits rebuilt from compact entries
		allBinsType() : _table((NUM_BINS + CHUNK_BINS - 1) / CHUNK_BINS, (Chunk *)NULL) {}
		~allBinsType() {
			for (size_t i = 0; i < _table.size(); i++) delete _table[i];
		}
		void stage(binNumType binNum, const binEntry &entry) {
			_stage.push_back(Staged(binNum, entry));
		}
		void finalize() {
			size_t n = _stage.size();
			for (size_t i = 0; i < n; i++) {
				binNumType b = _stage[i].bin;
				Chunk *&chunk = _table[b >> CHUNK_SHIFT];
				if (chunk == NULL) chunk = new Chunk();
				chunk->off[b & (CHUNK_BINS - 1)]++;
			}
			size_t base = 0;
			for (size_t c = 0; c < _table.size(); c++) {
				Chunk *chunk = _table[c];
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
			_entries.resize(n);
			for (size_t i = 0; i < n; i++) {
				binNumType b = _stage[i].bin;
				Chunk *chunk = _table[b >> CHUNK_SHIFT];
				_entries[chunk->base + chunk->off[b & (CHUNK_BINS - 1)]++] = _stage[i].entry;
			}
			vector<Staged>().swap(_stage);
			// one bit per bin: its entries are in nondecreasing start order
			for (size_t c = 0; c < _table.size(); c++) {
				Chunk *chunk = _table[c];
				if (chunk == NULL) continue;
				for (int b = 0; b < CHUNK_BINS; b++) {
					uint32_t lo = b == 0 ? 0 : chunk->off[b - 1];
					uint32_t hi = chunk->off[b];
					const binEntry *e = &_entries[chunk->base + lo];
					bool ok = true;
					for (uint32_t k = 1; k < hi - lo; k++) {
						if (e[k].start < e[k - 1].start) { ok = false; break; }
					}
					if (ok) chunk->sorted[b >> 6] |= (uint64_t)1 << (b & 63);
				}
			}
		}
		// false when the bin was never filled or is out of range
		bool find(binNumType binNum, const binEntry *&first, const binEntry *&last, bool &sorted) const {
			if ((uint64_t)binNum >= (uint64_t)NUM_BINS) return false;
			const Chunk *chunk = _table[binNum >> CHUNK_SHIFT];
			if (chunk == NULL) return false;
			int b = (int)(binNum & (CHUNK_BINS - 1));
			// after finalize(), off[b] is the end of bin b
			uint32_t lo = b == 0 ? 0 : chunk->off[b - 1];
			uint32_t hi = chunk->off[b];
			if (lo == hi) return false;
			first = &_entries[chunk->base + lo];
			last = first + (hi - lo);
			sorted = (chunk->sorted[b >> 6] >> (b & 63)) & 1;
			return true;
		}
	private:
		allBinsType(const allBinsType &);
		allBinsType &operator=(const allBinsType &);
		static const int CHUNK_SHIFT = 9;
		static const binNumType CHUNK_BINS = (binNumType)1 << CHUNK_SHIFT;
		struct Chunk {
			Chunk() : base(0) { memset(off, 0, sizeof(off)); memset(sorted, 0, sizeof(sorted)); }
			size_t base;
			uint64_t sorted[CHUNK_BINS / 64];
			uint32_t off[CHUNK_BINS];
		};
		struct Staged {
			Staged(binNumType b, const binEntry &e) : bin(b), entry(e) {}
			binNumType bin;
			binEntry entry;
		};
		vector<Chunk *> _table;
		vector<Staged> _stage;
		vector<binEntry> _entries;
	};
	typedef map<string, allBinsType *> mainMapType; //for each chrom, a map of bin num to RecordLists.
	mainMapType _mainMap;

	// last chrom looked up (queries) / inserted (load), to skip the map find.
	ChrRef _lastQueryChr;
	allBinsType *_lastQueryBins;
	bool _haveQueryChr;
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
	bool textEquals(int64_t off, const char *text, size_t len) const;

	bool addRecordToTree(Record *);
	bool stageCoords(binNumType startPos, binNumType endPos, int64_t meta);
	binNumType getBin(binNumType start, binNumType end) const;
	binNumType getBin(const Record *record) const;


};


#endif /* BINTREE_H_ */
