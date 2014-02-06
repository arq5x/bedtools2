/*
 * BinTree.h
 *
 *  Created on: Jan 5, 2013
 *      Author: nek3d
 */

#ifndef BINTREE_H_
#define BINTREE_H_

using namespace std;

#include <stdint.h>
#include <string>
#include <set>
#include <map>

#include "QuickString.h"
#include "RecordKeyList.h"
#include "ContextIntersect.h"

class FileRecordMgr;
class Record;

class BinTree {
public:
	BinTree(ContextIntersect *context);

	~BinTree();
	void loadDB();
	void getHits(Record *record, RecordKeyList &hitSet);

private:

	FileRecordMgr *_databaseFile;
	ContextIntersect *_context;

    //
    // BIN HANDLING
    //
	static const uint32_t NUM_BINS = 37450;
	static const uint32_t NUM_BIN_LEVELS = 7;

	// bins range in size from 16kb to 512Mb
	// Bin  0          spans 512Mbp,   # Level 1
	// Bins 1-8        span 64Mbp,     # Level 2
	// Bins 9-72       span 8Mbp,      # Level 3
	// Bins 73-584     span 1Mbp       # Level 4
	// Bins 585-4680   span 128Kbp     # Level 5
	// Bins 4681-37449 span 16Kbp      # Level 6
	uint32_t *_binOffsetsExtended;
	static const uint32_t _binFirstShift = 14;       /* How much to shift to get to finest bin. */
	static const uint32_t _binNextShift  = 3;        /* How much to shift to get to next larger bin. */

	typedef BTlist<const Record *> innerListType;
	typedef const BTlistNode<const Record *> * innerListIterType;
	typedef innerListType * binType;
	typedef binType * allBinsType;
	typedef QuickString mainKeyType;
	typedef map<mainKeyType, allBinsType> mainMapType;
	mainMapType _mainMap;

	bool _showBinMetrics;
	uint32_t _maxBinNumFound;
	map<uint32_t, int> _binsHit;

	bool addRecordToTree(const Record *);
	uint32_t getBin(uint32_t start, uint32_t end) const;
	uint32_t getBin(const Record *record) const;


};


#endif /* BINTREE_H_ */
