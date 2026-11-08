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


	typedef vector<Record *> binType;
	// For each bin number, a RecordList. Two-level paged table: a pointer per
	// CHUNK_BINS consecutive bin numbers, chunks allocated on first insert.
	class allBinsType {
	public:
		allBinsType() : _table((NUM_BINS + CHUNK_BINS - 1) / CHUNK_BINS, (binType *)NULL) {}
		~allBinsType() {
			for (size_t i = 0; i < _table.size(); i++) delete [] _table[i];
		}
		binType &operator[](binNumType binNum) {
			binType *&chunk = _table[binNum >> CHUNK_SHIFT];
			if (chunk == NULL) chunk = new binType[CHUNK_BINS];
			return chunk[binNum & (CHUNK_BINS - 1)];
		}
		// NULL when the bin was never filled or is out of range
		const binType *find(binNumType binNum) const {
			if ((uint64_t)binNum >= (uint64_t)NUM_BINS) return NULL;
			const binType *chunk = _table[binNum >> CHUNK_SHIFT];
			return chunk == NULL ? NULL : &chunk[binNum & (CHUNK_BINS - 1)];
		}
	private:
		allBinsType(const allBinsType &);
		allBinsType &operator=(const allBinsType &);
		static const int CHUNK_SHIFT = 9;
		static const binNumType CHUNK_BINS = (binNumType)1 << CHUNK_SHIFT;
		vector<binType *> _table;
	};
	typedef map<string, allBinsType *> mainMapType; //for each chrom, a map of bin num to RecordLists.
	mainMapType _mainMap;

	map<binNumType, int> _binsHit;

	bool addRecordToTree(Record *);
	binNumType getBin(binNumType start, binNumType end) const;
	binNumType getBin(const Record *record) const;


};


#endif /* BINTREE_H_ */
