#include "BinTree.h"
#include "FileRecordMgr.h"


BinTree::BinTree(ContextIntersect *context) 
:  _context(context), _lastQueryBins(NULL), _haveQueryChr(false), _lastAddBins(NULL) {
}

BinTree::~BinTree() {
	for (mainMapType::iterator it = _mainMap.begin(); it != _mainMap.end(); ++it) {
		delete it->second;
	}
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
	if (!_haveQueryChr || chr != _lastQueryChr) {
		mainMapType::iterator mainIter = _mainMap.find(chr);
		_lastQueryChr = chr;
		_haveQueryChr = true;
		_lastQueryBins = (mainIter == _mainMap.end()) ? NULL : mainIter->second;
	}
	if (_lastQueryBins == NULL) {
		//given chrom not even in map.
		return;
	}

    binNumType startPos = record->getStartPos();
    binNumType endPos = record->getEndPos();

    binNumType startBin = (startPos >> _binFirstShift);
    binNumType endBin = ((endPos-1) >> _binFirstShift);


	const allBinsType &bins = *_lastQueryBins;

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
        	const binType *binPtr = bins.find(j);
        	if (binPtr == NULL) {
        		continue;
        	}
        	const binType &bin = *binPtr;

        	for (binType::const_iterator iter = bin.begin(); iter != bin.end(); iter++) {
            	// Record::sameChromIntersects returns false when minEnd < maxStart.
            	if (iter->end < startPos || iter->start > endPos) {
            		continue;
            	}
            	Record *dbRec = iter->record;
            	if (record->intersects(dbRec,
                                       _context->getSameStrand(),
                                       _context->getDiffStrand(),
            			               _context->getOverlapFractionA(),
                                       _context->getOverlapFractionB(),
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
		fprintf(stderr, "ERROR: Received illegal bin number %" PRId_BINNUMTYPE " from getBin call.\n"
                                "Maximum values is: %" PRId_BINNUMTYPE "\n"
                                "This typically means that your coordinates are\n"
                                "negative or too large to represent in the data\n"
                                "structure bedtools uses to find intersections.", binNum, NUM_BINS);
		return false;
	}
	if (_lastAddBins == NULL || chr != _lastAddChr) {
		allBinsType *&slot = _mainMap[chr];
		if (slot == NULL) {
			slot = new allBinsType();
		}
		_lastAddChr = chr;
		_lastAddBins = slot;
	}
	allBinsType *bins = _lastAddBins;
	_haveQueryChr = false;
	binEntry entry = {startPos, endPos, record};
	(*bins)[binNum].push_back(entry);
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
