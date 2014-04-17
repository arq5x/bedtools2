#include "BinTree.h"
#include "FileRecordMgr.h"


BinTree::BinTree(ContextIntersect *context)
: _databaseFile(NULL),
  _context(context),
  _binOffsetsExtended(NULL),
  _showBinMetrics(false),
  _maxBinNumFound(0)
 {
	_binOffsetsExtended = new uint32_t[NUM_BIN_LEVELS];
	memset(_binOffsetsExtended, 0, NUM_BIN_LEVELS * sizeof(uint32_t));

	//start at idx 1, because the memset above already initialized
	//the first idx to zero, which is what we want.
	for (uint32_t i= 1; i < NUM_BIN_LEVELS; i++) {
		_binOffsetsExtended[i] = _binOffsetsExtended[i-1] + (1 << ((NUM_BIN_LEVELS - i -1) * 3));
	}
}

BinTree::~BinTree() {
	//TBD: pass all elements in tree back to FRM for proper cleanup/deletion
	for (mainMapType::iterator mainIter = _mainMap.begin(); mainIter != _mainMap.end(); mainIter++) {
		allBinsType bins = mainIter->second;
		if (bins == NULL) {
			fprintf(stderr, "ERROR: In BinTree destructor: found chromosome with NULL bin array.\n");
			continue;
		}
		if (!_showBinMetrics) { //don't clean up bins when simply reporting metrics.

			for (uint32_t i=0; i < NUM_BINS; i++) {
				binType bin = bins[i];
				if (bin == NULL) {
					continue;
				}
				for (innerListIterType listIter = bin->begin(); listIter != bin->end(); listIter = bin->next()) {
					const Record *record = listIter->value();
					_databaseFile->deleteRecord(record);
				}
				delete bin;
				bin = NULL;
			}
		}
		delete [] bins;
		bins = NULL;
	}
	delete [] _binOffsetsExtended;

	if (_showBinMetrics) {
		map<int, int> hitsHistogram;
		FILE *fp = fopen("BinsHitFile.txt", "w");
		fprintf(fp, "The largest bin was %u\n", _maxBinNumFound);
		fprintf(fp, "There were %d different bins hit by database.\n", (int)_binsHit.size());
		for (map<uint32_t, int>::iterator binIter = _binsHit.begin(); binIter != _binsHit.end(); binIter++) {
			uint32_t binNum = binIter->first;
			int numHits = binIter->second;
			fprintf(fp, "%u\t%d\n", binNum, numHits);
			hitsHistogram[numHits]++;
		}
		fclose(fp);
		fp = fopen("BinHitsHistogram.txt", "w");
		fprintf(fp, "NumHits\tNumBins\n");
		for (map<int, int>::iterator histIter = hitsHistogram.begin(); histIter != hitsHistogram.end(); histIter++) {
			fprintf(fp, "%d\t%d\n", histIter->first, histIter->second);
		}
		fclose(fp);
	}
}

void BinTree::loadDB()
{
	_databaseFile = _context->getFile(_context->getDatabaseFileIdx());

	Record *record = NULL;
	while (!_databaseFile->eof()) {
		record = _databaseFile->getNextRecord();
		//In addition to NULL records, we also don't want to add unmapped reads.
		if (record == NULL || record->isUnmapped()) {
			continue;
		}

		if (!addRecordToTree(record)) {
			fprintf(stderr, "ERROR: Unable to add record to tree.\n");
			_databaseFile->close();
			exit(1);
		}
	}
}

void BinTree::getHits(Record *record, RecordKeyList &hitSet)
{
	if (_showBinMetrics) {
		return; //don't care about query entries just yet.
	}
	if (record->isUnmapped()) {
		return;
	}
    const QuickString &chr = record->getChrName();
	mainMapType::iterator mainIter = _mainMap.find(chr);
	if (mainIter == _mainMap.end()) {
		//given chrom not even in map.
		return;
	}

    uint32_t startBin = 0;
    uint32_t endBin = 0;

    uint32_t startPos = record->getStartPos();
    uint32_t endPos = record->getEndPos();

    startBin = (startPos >> _binFirstShift);
    endBin = ((endPos-1) >> _binFirstShift);


	const allBinsType bins = mainIter->second;

    /* SYNOPSIS:
         1. We loop through each UCSC BIN level for feature A's chrom.
         2. For each BIN, we loop through each B feature and add it to
            hits if it meets all of the user's requests, which include:
               (a) overlap fractio, (b) strandedness, (c) reciprocal overlap
    */
    for (uint32_t i = 0; i < NUM_BIN_LEVELS; i++) {
        uint32_t offset = _binOffsetsExtended[i];
        for (uint32_t j = (startBin+offset); j <= (endBin+offset); j++)  {

        	// move to the next bin if this one is empty
        	const binType &bin = bins[j];
        	if (bin == NULL) {
        		//no list of records in this bin.
        		continue;
        	}
        	if (bin->empty()) {
        		//bin has list, but it's empty.
        		//Actually, this should never happen, so throw an error.
        		fprintf(stderr, "ERROR: found empty list for bin %u of chromosome %s\n",
        				j, chr.c_str());
        		continue;
        	}
            for (innerListIterType listIter = bin->begin(); listIter != bin->end(); listIter = bin->next()) {
            	const Record *dbRec = listIter->value();
            	if (record->intersects(dbRec, _context->getSameStrand(), _context->getDiffStrand(),
            			_context->getOverlapFraction(), _context->getReciprocal())) {
            		hitSet.push_back(dbRec);
            	}
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}

bool BinTree::addRecordToTree(const Record *record)
{
	// Get chr, bin. allocate all bins and single bins as needed.
	const QuickString &chr = record->getChrName();
	uint32_t startPos = (uint32_t)(record->getStartPos());
	uint32_t endPos = (uint32_t)(record->getEndPos());

	//is this chr currently in the main map?
	allBinsType bins = NULL;
	mainMapType::iterator mainIter = _mainMap.find(chr);
	if (mainIter == _mainMap.end()) {
		//this is a new chr NOT currently in the map.
		bins = new binType[NUM_BINS];
		memset(bins, 0, NUM_BINS * sizeof(binType));
		_mainMap[chr] = bins;
	} else {
		bins = mainIter->second;
	}
	uint32_t binNum = getBin(startPos, endPos);

	if (_showBinMetrics) {
		if (binNum > _maxBinNumFound) {
			_maxBinNumFound = binNum;
		}
		_binsHit[binNum]++;
		return true;
	}

	if (binNum < 0 || binNum >= NUM_BINS) {
		fprintf(stderr, "ERROR: Received illegal bin number %u from getBin call.\n", binNum);
		return false;
	}
	binType &bin = bins[binNum];
	if (bin == NULL) {
		bin = new innerListType;
	}
	bin->push_back(record);

	return true;
}

uint32_t BinTree::getBin(const Record *record) const {
	return getBin((uint32_t)(record->getStartPos()), (uint32_t)(record->getEndPos()));
}

uint32_t BinTree::getBin(uint32_t start, uint32_t end) const {
    --end;
    start >>= _binFirstShift;
    end   >>= _binFirstShift;

    for (uint32_t i = 0; i < NUM_BIN_LEVELS; ++i) {
        if (start == end) {
        	return _binOffsetsExtended[i] + start;
        }
        start >>= _binNextShift;
        end   >>= _binNextShift;
    }
    //failure
    return -1;
}
