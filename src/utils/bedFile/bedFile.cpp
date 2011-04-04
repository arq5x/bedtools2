/*****************************************************************************
  bedFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"


/************************************************
Helper functions
*************************************************/
void splitBedIntoBlocks(const BED &bed, int lineNum, bedVector &bedBlocks) {

    if (bed.otherFields.size() < 6) {
        cerr << "Input error: Cannot split into blocks. Found interval with fewer than 12 columns on line " << lineNum << "." << endl;
        exit(1);
    }

    int blockCount = atoi(bed.otherFields[3].c_str());
    if ( blockCount <= 0 ) {
        cerr << "Input error: found interval having <= 0 blocks on line " << lineNum << "." << endl;
        exit(1);
    }
    else if ( blockCount == 1 ) {
        //take a short-cut for single blocks
        bedBlocks.push_back(bed);
    }
    else {
        // get the comma-delimited strings for the BED12 block starts and block ends.
        string blockSizes(bed.otherFields[4]);
        string blockStarts(bed.otherFields[5]);

        vector<int> sizes;
        vector<int> starts;
        Tokenize(blockSizes, sizes, ",");
        Tokenize(blockStarts, starts, ",");

        if ( sizes.size() != (size_t) blockCount || starts.size() != (size_t) blockCount ) {
            cerr << "Input error: found interval with block-counts not matching starts/sizes on line " << lineNum << "." << endl;
            exit(1);
        }

        // add each BED block to the bedBlocks vector
        for (UINT i = 0; i < (UINT) blockCount; ++i) {
            CHRPOS blockStart = bed.start + starts[i];
            CHRPOS blockEnd   = bed.start + starts[i] + sizes[i];
            BED currBedBlock(bed.chrom, blockStart, blockEnd, bed.name, bed.score, bed.strand, bed.otherFields);
            bedBlocks.push_back(currBedBlock);
        }
    }
}


/***********************************************
Sorting comparison functions
************************************************/
bool sortByChrom(BED const &a, BED const &b) {
    if (a.chrom < b.chrom) return true;
    else return false;
};

bool sortByStart(const BED &a, const BED &b) {
    if (a.start < b.start) return true;
    else return false;
};

bool sortBySizeAsc(const BED &a, const BED &b) {

    CHRPOS aLen = a.end - a.start;
    CHRPOS bLen = b.end - b.start;

    if (aLen < bLen) return true;
    else return false;
};

bool sortBySizeDesc(const BED &a, const BED &b) {

    CHRPOS aLen = a.end - a.start;
    CHRPOS bLen = b.end - b.start;

    if (aLen > bLen) return true;
    else return false;
};

bool sortByScoreAsc(const BED &a, const BED &b) {
    if (a.score < b.score) return true;
    else return false;
};

bool sortByScoreDesc(const BED &a, const BED &b) {
    if (a.score > b.score) return true;
    else return false;
};

bool byChromThenStart(BED const &a, BED const &b) {

    if (a.chrom < b.chrom) return true;
    else if (a.chrom > b.chrom) return false;

    if (a.start < b.start) return true;
    else if (a.start >= b.start) return false;

    return false;
};


/*******************************************
Class methods
*******************************************/

// Constructor
BedFile::BedFile(string &bedFile)
: bedFile(bedFile),
  _typeIsKnown(false)
{}

// Destructor
BedFile::~BedFile(void) {
}


void BedFile::Open(void) {
    if (bedFile == "stdin") {
        _bedStream = &cin;
    }
    // New method thanks to Assaf Gordon
    else if ((isGzipFile(bedFile) == false) && (isRegularFile(bedFile) == true)) {
       // open an ifstream
        ifstream beds(bedFile.c_str(), ios::in);

        // can we open the file?
        if ( !beds ) {
            cerr << "Error: The requested bed file (" << bedFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        else {
            // if so, close it (this was just a test)
            beds.close();
            // now set a pointer to the stream so that we
            _bedStream = new ifstream(bedFile.c_str(), ios::in);
        }
    }
    else if ((isGzipFile(bedFile) == true) && (isRegularFile(bedFile) == true)) {
        igzstream beds(bedFile.c_str(), ios::in);
        if ( !beds ) {
            cerr << "Error: The requested bed file (" << bedFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        else {
            // if so, close it (this was just a test)
            beds.close();
            // now set a pointer to the stream so that we
            _bedStream = new igzstream(bedFile.c_str(), ios::in);
        }
    }
    else {
        cerr << "Error: Unexpected file type (" << bedFile << "). Exiting!" << endl;
        exit(1);
    }
}


// Close the BED file
void BedFile::Close(void) {
    if (bedFile != "stdin") delete _bedStream;
}


BedLineStatus BedFile::GetNextBed(BED &bed, int &lineNum) {

    // make sure there are still lines to process.
    // if so, tokenize, validate and return the BED entry.
    if (_bedStream->good()) {
        string bedLine;
        vector<string> bedFields;
        bedFields.reserve(12);

        // parse the bedStream pointer
        getline(*_bedStream, bedLine);
        lineNum++;

        // split into a string vector.
        Tokenize(bedLine,bedFields);

        // load the BED struct as long as it's a valid BED entry.
        return parseLine(bed, bedFields, lineNum);
    }

    // default if file is closed or EOF
    return BED_INVALID;
}


void BedFile::FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end,
                                 string strand, vector<BED> &hits, bool forceStrand) {

    BIN startBin, endBin;
    startBin = (start >> _binFirstShift);
    endBin = ((end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                if (overlaps(bedItr->start, bedItr->end, start, end) > 0) {
                    // skip the hit if not on the same strand (and we care)
                    if (forceStrand == false) hits.push_back(*bedItr);
                    else if ( (forceStrand == true) && (strand == bedItr->strand)) {
                         hits.push_back(*bedItr);
                    }
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}


bool BedFile::FindOneOrMoreOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand,
    bool forceStrand, float overlapFraction) {

    BIN startBin, endBin;
    startBin = (start   >> _binFirstShift);
    endBin   = ((end-1) >> _binFirstShift);

    CHRPOS aLength = (end - start);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[chrom][j].end();
            for (; bedItr != bedEnd; ++bedItr) {

                CHRPOS s = max(start, bedItr->start);
                CHRPOS e = min(end, bedItr->end);
                // the number of overlapping bases b/w a and b
                int overlapBases = (e - s);

                // do we have sufficient overlap?
                if ( (float) overlapBases / (float) aLength  >= overlapFraction) {
                    // skip the hit if not on the same strand (and we care)
                    if (forceStrand == false) return true;
                    else if ( (forceStrand == true) && (strand == bedItr->strand)) {
                        return true;
                    }
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return false;
}


bool BedFile::FindOneOrMoreReciprocalOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand,
    bool forceStrand, float overlapFraction) {

    BIN startBin, endBin;
    startBin = (start >> _binFirstShift);
    endBin = ((end-1) >> _binFirstShift);

    CHRPOS aLength = (end - start);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[chrom][j].end();
            for (; bedItr != bedEnd; ++bedItr) {
                CHRPOS s = max(start, bedItr->start);
                CHRPOS e = min(end, bedItr->end);

                // the number of overlapping bases b/w a and b
                int overlapBases = (e - s);

                // do we have sufficient overlap?
                if ( (float) overlapBases / (float) aLength  >= overlapFraction) {
                    CHRPOS bLength = (bedItr->end - bedItr->start);
                    float bOverlap = ( (float) overlapBases / (float) bLength );
                    if ((forceStrand == false) && (bOverlap >= overlapFraction)) {
                        return true;
                    }
                    else if ( (forceStrand == true) && (strand == bedItr->strand) && (bOverlap >= overlapFraction)) {
                        return true;
                    }
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return false;
}


void BedFile::countHits(const BED &a, bool forceStrand) {

    BIN startBin, endBin;
    startBin = (a.start >> _binFirstShift);
    endBin = ((a.end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j) {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to
            // the list of hits.
            vector<BEDCOV>::iterator bedItr = bedCovMap[a.chrom][j].begin();
            vector<BEDCOV>::iterator bedEnd = bedCovMap[a.chrom][j].end();
            for (; bedItr != bedEnd; ++bedItr) {

                // skip the hit if not on the same strand (and we care)
                if (forceStrand && (a.strand != bedItr->strand)) {
                    continue;
                }
                else if (overlaps(bedItr->start, bedItr->end, a.start, a.end) > 0) {

                    bedItr->count++;
                    if (a.zeroLength == false) {
                        bedItr->depthMap[a.start+1].starts++;
                        bedItr->depthMap[a.end].ends++;
                    }
                    else {
                        // correct for the fact that we artificially expanded the zeroLength feature
                        bedItr->depthMap[a.start+2].starts++;
                        bedItr->depthMap[a.end-1].ends++;                        
                    }

                    if (a.start < bedItr->minOverlapStart) {
                        bedItr->minOverlapStart = a.start;
                    }
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}


void BedFile::countSplitHits(const vector<BED> &bedBlocks, bool forceStrand) {

    // set to track the distinct B features that had coverage.
    // we'll update the counts of coverage for these features by one
    // at the end of this function to avoid over-counting.
    set< vector<BEDCOV>::iterator > validHits;

    vector<BED>::const_iterator blockItr  = bedBlocks.begin();
    vector<BED>::const_iterator blockEnd  = bedBlocks.end();
    for (; blockItr != blockEnd; ++blockItr) {

        BIN startBin, endBin;
        startBin = (blockItr->start >> _binFirstShift);
        endBin = ((blockItr->end-1) >> _binFirstShift);

        // loop through each bin "level" in the binning hierarchy
        for (BINLEVEL i = 0; i < _binLevels; ++i) {

            // loop through each bin at this level of the hierarchy
            BIN offset = _binOffsetsExtended[i];
            for (BIN j = (startBin+offset); j <= (endBin+offset); ++j) {

                // loop through each feature in this chrom/bin and see if it overlaps
                // with the feature that was passed in.  if so, add the feature to
                // the list of hits.
                vector<BEDCOV>::iterator bedItr = bedCovMap[blockItr->chrom][j].begin();
                vector<BEDCOV>::iterator bedEnd = bedCovMap[blockItr->chrom][j].end();
                for (; bedItr != bedEnd; ++bedItr) {

                    // skip the hit if not on the same strand (and we care)
                    if (forceStrand && (blockItr->strand != bedItr->strand)) {
                        continue;
                    }
                    else if (overlaps(bedItr->start, bedItr->end, blockItr->start, blockItr->end) > 0) {
                        if (blockItr->zeroLength == false) {
                            bedItr->depthMap[blockItr->start+1].starts++;
                            bedItr->depthMap[blockItr->end].ends++;
                        }
                        else {
                            // correct for the fact that we artificially expanded the zeroLength feature
                            bedItr->depthMap[blockItr->start+2].starts++;
                            bedItr->depthMap[blockItr->end-1].ends++;
                        }

                        validHits.insert(bedItr);
                        if (blockItr->start < bedItr->minOverlapStart)
                            bedItr->minOverlapStart = blockItr->start;
                    }
                }
            }
            startBin >>= _binNextShift;
            endBin >>= _binNextShift;
        }
    }
    // incrment the count of overlap by one for each B feature that overlapped
    // the current passed hit.  This is necessary to prevent over-counting for
    // each "split"" of a single read.
    set< vector<BEDCOV>::iterator >::iterator validHitsItr = validHits.begin();
    set< vector<BEDCOV>::iterator >::iterator validHitsEnd = validHits.end();
    for (; validHitsItr != validHitsEnd; ++validHitsItr)
        // the validHitsItr points to another itr, hence the (*itr)-> dereferencing.
        // ugly, but that's C++.
        (*validHitsItr)->count++;
}


void BedFile::countListHits(const BED &a, int index, bool forceStrand) {

    BIN startBin, endBin;
    startBin = (a.start >> _binFirstShift);
    endBin = ((a.end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j) {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to
            // the list of hits.
            vector<BEDCOVLIST>::iterator bedItr = bedCovListMap[a.chrom][j].begin();
            vector<BEDCOVLIST>::iterator bedEnd = bedCovListMap[a.chrom][j].end();
            for (; bedItr != bedEnd; ++bedItr) {

                if (forceStrand && (a.strand != bedItr->strand)) {
                    continue;
                }
                else if (overlaps(bedItr->start, bedItr->end, a.start, a.end) > 0) {
                    bedItr->counts[index]++;
                    if (a.zeroLength == false) {
                        bedItr->depthMapList[index][a.start+1].starts++;
                        bedItr->depthMapList[index][a.end].ends++;
                    }
                    else {
                        // correct for the fact that we artificially expanded the zeroLength feature
                        bedItr->depthMapList[index][a.start+2].starts++;
                        bedItr->depthMapList[index][a.end-1].ends++;                        
                    }

                    if (a.start < bedItr->minOverlapStarts[index]) {
                        bedItr->minOverlapStarts[index] = a.start;
                    }
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}

void BedFile::setZeroBased(bool zeroBased) { this->isZeroBased = zeroBased; }

void BedFile::setGff (bool gff) { this->_isGff = gff; }


void BedFile::setVcf (bool vcf) { this->_isVcf = vcf; }


void BedFile::setFileType (FileType type) {
    _fileType    = type;
    _typeIsKnown = true;
}


void BedFile::setBedType (int colNums) {
    bedType = colNums;
}


void BedFile::loadBedFileIntoMap() {

    BED bedEntry, nullBed;
    int lineNum = 0;
    BedLineStatus bedStatus;

    Open();
    while ((bedStatus = GetNextBed(bedEntry, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            BIN bin = getBin(bedEntry.start, bedEntry.end);
            bedMap[bedEntry.chrom][bin].push_back(bedEntry);
            bedEntry = nullBed;
        }
    }
    Close();
}


void BedFile::loadBedCovFileIntoMap() {

    BED bedEntry, nullBed;
    int lineNum = 0;
    BedLineStatus bedStatus;

    Open();
    while ((bedStatus = GetNextBed(bedEntry, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            BIN bin = getBin(bedEntry.start, bedEntry.end);

            BEDCOV bedCov;
            bedCov.chrom        = bedEntry.chrom;
            bedCov.start        = bedEntry.start;
            bedCov.end          = bedEntry.end;
            bedCov.name         = bedEntry.name;
            bedCov.score        = bedEntry.score;
            bedCov.strand       = bedEntry.strand;
            bedCov.otherFields  = bedEntry.otherFields;
            bedCov.zeroLength   = bedEntry.zeroLength;
            bedCov.count = 0;
            bedCov.minOverlapStart = INT_MAX;

            bedCovMap[bedEntry.chrom][bin].push_back(bedCov);
            bedEntry = nullBed;
        }
    }
    Close();
}

void BedFile::loadBedCovListFileIntoMap() {

    BED bedEntry, nullBed;
    int lineNum = 0;
    BedLineStatus bedStatus;

    Open();
    while ((bedStatus = GetNextBed(bedEntry, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            BIN bin = getBin(bedEntry.start, bedEntry.end);

            BEDCOVLIST bedCovList;
            bedCovList.chrom        = bedEntry.chrom;
            bedCovList.start        = bedEntry.start;
            bedCovList.end          = bedEntry.end;
            bedCovList.name         = bedEntry.name;
            bedCovList.score        = bedEntry.score;
            bedCovList.strand       = bedEntry.strand;
            bedCovList.otherFields  = bedEntry.otherFields;
            bedCovList.zeroLength   = bedEntry.zeroLength;

            bedCovListMap[bedEntry.chrom][bin].push_back(bedCovList);
            bedEntry = nullBed;
        }
    }
    Close();
}


void BedFile::loadBedFileIntoMapNoBin() {

    BED bedEntry, nullBed;
    int lineNum = 0;
    BedLineStatus bedStatus;

    Open();
    while ((bedStatus = this->GetNextBed(bedEntry, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            bedMapNoBin[bedEntry.chrom].push_back(bedEntry);
            bedEntry = nullBed;
        }
    }
    Close();

    // sort the BED entries for each chromosome
    // in ascending order of start position
    for (masterBedMapNoBin::iterator m = this->bedMapNoBin.begin(); m != this->bedMapNoBin.end(); ++m) {
        sort(m->second.begin(), m->second.end(), sortByStart);
    }
}
