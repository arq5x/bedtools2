/*
 * BinTree.h
 *
 *  Created on: Jan 5, 2013
 *      Author: nek3d
 */

#ifndef BINTREE_H_
#define BINTREE_H_

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
	typedef int32_t binNumType;
	static const binNumType NUM_BINS = 37450;
	static const binNumType NUM_BIN_LEVELS = 7;

	// bins range in size from 16kb to 512Mb
	// Bin  0          spans 512Mbp,   # Level 1
	// Bins 1-8        span 64Mbp,     # Level 2
	// Bins 9-72       span 8Mbp,      # Level 3
	// Bins 73-584     span 1Mbp       # Level 4
	// Bins 585-4680   span 128Kbp     # Level 5
	// Bins 4681-37449 span 16Kbp      # Level 6
	binNumType *_binOffsetsExtended;
	static const binNumType _binFirstShift = 14;       /* How much to shift to get to finest bin. */
	static const binNumType _binNextShift  = 3;        /* How much to shift to get to next larger bin. */

//	typedef RecordList innerListType;
//	typedef const RecordListNode * innerListIterType;
//	typedef innerListType * binType;
//	typedef binType * allBinsType;
//	typedef string mainKeyType;
//	typedef map<mainKeyType, allBinsType> mainMapType;
//	mainMapType _mainMap;

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
