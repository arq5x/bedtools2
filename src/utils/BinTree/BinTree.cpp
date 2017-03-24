#include "BinTree.h"
#include "FileRecordMgr.h"


BinTree::BinTree(ContextIntersect *context)
:  _context(context),
  _binOffsetsExtended(NULL)
 {
	_binOffsetsExtended = new binNumType[NUM_BIN_LEVELS];
	memset(_binOffsetsExtended, 0, NUM_BIN_LEVELS * sizeof(binNumType));

	//start at idx 1, because the memset above already initialized
	//the first idx to zero, which is what we want.
	for (binNumType i= 1; i < NUM_BIN_LEVELS; i++) {
		_binOffsetsExtended[i] = _binOffsetsExtended[i-1] + (1 << ((NUM_BIN_LEVELS - i -1) * 3));
	}
}

BinTree::~BinTree() {
	delete [] _binOffsetsExtended;
}

void BinTree::loadDB()
{
	for (int i=0; i < _context->getNumDatabaseFiles(); i++) {
		FileRecordMgr *databaseFile = _context->getDatabaseFile(i);

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
}

void BinTree::getHits(Record *record, RecordKeyVector &hitSet)
{
	if (record->isUnmapped()) {
		return;
	}
    const string &chr = record->getChrName();
	mainMapType::iterator mainIter = _mainMap.find(chr);
	if (mainIter == _mainMap.end()) {
		//given chrom not even in map.
		return;
	}

    binNumType startPos = record->getStartPos();
    binNumType endPos = record->getEndPos();

    binNumType startBin = (startPos >> _binFirstShift);
    binNumType endBin = ((endPos-1) >> _binFirstShift);


	allBinsType &bins = mainIter->second;

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
        	allBinsType::iterator allBinsIter = bins.find(j);
        	if (allBinsIter == bins.end()) {
        		continue;
        	}
        	binType &bin = allBinsIter->second;

        	for (binType::iterator iter = bin.begin(); iter != bin.end(); iter++) {
            	Record *dbRec = *iter;
            	if (record->intersects(dbRec,
                                       _context->getSameStrand(),
                                       _context->getDiffStrand(),
            			               _context->getOverlapFractionA(),
                                       _context->getOverlapFractionB(),
                                       _context->getReciprocalFraction(),
                                       _context->getEitherFraction()
                                      )
                    )
                {
            		hitSet.push_back(dbRec);
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

bool BinTree::addRecordToTree(Record *record)
{
	// Get chr, bin.
	const string &chr = record->getChrName();
	binNumType startPos = (binNumType)(record->getStartPos());
	binNumType endPos = (binNumType)(record->getEndPos());
	binNumType binNum = getBin(startPos, endPos);

	if (binNum < 0 || binNum >= NUM_BINS) {
		fprintf(stderr, "ERROR: Received illegal bin number %u from getBin call.\n", binNum);
		return false;
	}
	_mainMap[chr][binNum].push_back(record);
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
