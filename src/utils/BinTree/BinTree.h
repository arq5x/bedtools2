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
	typedef map<binNumType, binType> allBinsType; //for each bin number, have a RecordList
	typedef map<string, allBinsType> mainMapType; //for each chrom, a map of bin num to RecordLists.
	mainMapType _mainMap;

	map<binNumType, int> _binsHit;

	bool addRecordToTree(Record *);
	binNumType getBin(binNumType start, binNumType end) const;
	binNumType getBin(const Record *record) const;


};


#endif /* BINTREE_H_ */
