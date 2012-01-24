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
  _isGff(false),
  _isVcf(false),
  _typeIsKnown(false),
  _merged_start(-1),
  _merged_end(-1),
  _merged_chrom(""),
  _prev_start(-1),
  _prev_chrom("")
{}

// Destructor
BedFile::~BedFile(void) {
}


void BedFile::Open(void) {
    
    _bedFields.reserve(12);
    
    if (bedFile == "stdin" || bedFile == "-") {
        _bedStream = &cin;
    }
    else {
        _bedStream = new ifstream(bedFile.c_str(), ios::in);
        
        if( isGzipFile(_bedStream) ) {
            delete _bedStream;
            _bedStream = new igzstream(bedFile.c_str(), ios::in);
        }
        if ( !(_bedStream->good()) ) {
            cerr << "Error: The requested bed file (" << bedFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
    }
    // save the file's header (if there is one)
    GetHeader();
}

// Rewind the pointer back to the beginning of the file
void BedFile::Rewind(void) {
    _bedStream->seekg(0, ios::beg);
}

// Jump to a specific byte in the file
void BedFile::Seek(unsigned long offset) {
    _bedStream->seekg(offset);
}

// Jump to a specific byte in the file
bool BedFile::Empty(void) {
    return _status == BED_INVALID || _status == BED_BLANK;
}

// Close the BED file
void BedFile::Close(void) {
    if (bedFile != "stdin" && bedFile != "-") delete _bedStream;
}

void BedFile::GetLine(void) {
    // parse the bedStream pointer
    getline(*_bedStream, _bedLine);
    // increment the line number
    _lineNum++;
    // split into a string vector.
    Tokenize(_bedLine, _bedFields);
}

// Extract and store the header for the file.
void BedFile::GetHeader(void) {
    while(getline(*_bedStream, _bedLine))
    {
        _lineNum++;
        // look for header lines.  ^# headers can span multiple lines, 
        // but ^[browser|track|chrom] headers must occur on the 1st line.
        if ( (_bedLine.find("#")       == 0) ||
             (_bedLine.find("browser") == 0) ||
             (_bedLine.find("track")   == 0) 
           )
        {
            _header += _bedLine + '\n';
        }
        // we are done with the header. stop looking
        // and indicate that the first data line has been read
        // (i.e., _bedLine now houses the first data line)
        else
        {
            _firstLine = true;
            break;
        }
    }
}

// Dump the header
void BedFile::PrintHeader(void) {
    cout << _header;
}


bool BedFile::GetNextBed(BED &bed, bool forceSorted) {

    // make sure there are still lines to process.
    // if so, tokenize, validate and return the BED entry.
    _bedFields.clear();
    // clear out the previous bed's data
    if (_bedStream->good()) {
        // read the next line in the file and parse into discrete fields
        if (!_firstLine)
            GetLine();
        else {
            // handle the first line as a special case because
            // of reading the header.
            Tokenize(_bedLine, _bedFields);
            _firstLine = false;
        }
        // load the BED struct as long as it's a valid BED entry.
        _status = parseLine(bed, _bedFields);
        if (_status == BED_INVALID) return false;
        
        if (_status == BED_VALID) {
            if (bed.chrom == _prev_chrom) {
                if ((int) bed.start >= _prev_start) {
                    _prev_chrom = bed.chrom;
                    _prev_start = bed.start;
                }
                else if (forceSorted) {
                    cerr << "ERROR: input file: (" << bedFile 
                         << ") is not sorted by chrom then start." << endl
                         << "       The start coordinate at line " << _lineNum 
                         << " is less than the start at line " << _lineNum-1 
                         << endl;
                    exit(1);
                }
            }
            else if (bed.chrom != _prev_chrom) {
                _prev_chrom = bed.chrom;
                _prev_start = bed.start;
            }
            return true;
        }
        else if (_status == BED_HEADER || _status == BED_BLANK) 
        {
            return true;
        }
    }

    // default if file is closed or EOF
    _status = BED_INVALID;
    return false;
}


bool BedFile::GetNextMergedBed(BED &merged_bed) {

    if (_bedStream->good()) {
        BED bed;
        // force sorting; hence third param = true
        while (GetNextBed(bed, true)) {
            if (_status == BED_VALID) {
                if (((int) bed.start - _merged_end > 0) || 
                   (_merged_end < 0) || 
                   (bed.chrom != _merged_chrom))
                {
                    if (_merged_start >= 0) {
                        merged_bed.chrom = _merged_chrom;
                        merged_bed.start = _merged_start;
                        merged_bed.end   = _merged_end;

                        _merged_chrom = bed.chrom;
                        _merged_start = bed.start;
                        _merged_end   = bed.end;

                        return true;
                    }
                    else {
                        _merged_start = bed.start;
                        _merged_chrom = bed.chrom;
                        _merged_end = bed.end;
                    }
                }
                else if ((int) bed.end > _merged_end)
                {   
                    _merged_end = bed.end;
                }
            }
        }
        // handle the last merged block in the file.
        if (_status == BED_INVALID)
        {
            merged_bed.chrom = _merged_chrom;
            merged_bed.start = _merged_start;
            merged_bed.end   = _merged_end;
            return true;
        }
    }
    return false;
}


void BedFile::allHits(string chrom, CHRPOS start, CHRPOS end, string strand, 
                                 vector<BED> &hits, bool sameStrand, bool diffStrand,
                                 float overlapFraction, bool reciprocal) {

    BIN startBin, endBin;
    startBin = (start >> _binFirstShift);
    endBin = ((end-1) >> _binFirstShift);
    CHRPOS aLength = (end - start);

    /* SYNOPSIS:
         1. We loop through each UCSC BIN level for feature A's chrom.
         2. For each BIN, we loop through each B feature and add it to
            hits if it meets all of the user's requests, which include:
               (a) overlap fractio, (b) strandedness, (c) reciprocal overlap
    */
    for (BINLEVEL i = 0; i < _binLevels; ++i) {
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {
            // move to the next bin if this one is empty
            if (bedMap[chrom][j].empty()) continue;
            vector<BED>::const_iterator bedItr = bedMap[chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[chrom][j].end();
            for (; bedItr != bedEnd; ++bedItr) {
                CHRPOS s = max(start, bedItr->start);
                CHRPOS e = min(end, bedItr->end);
                int overlapBases = (e - s); 
                // 1. is there sufficient overlap w.r.t A?
                if ( (float) overlapBases / (float) aLength  >= overlapFraction) {
                    CHRPOS bLength = (bedItr->end - bedItr->start);
                    float bOverlap = ( (float) overlapBases / (float) bLength );
                    bool strands_are_same = (strand == bedItr->strand);
                    // 2. does the overlap meet the user's strand repuirements?
                    if ( (sameStrand == false && diffStrand == false)
                         ||
                         (sameStrand == true && strands_are_same == true)
                         ||
                         (diffStrand == true && strands_are_same == false)
                       )
                    {
                        // 3. did the user request reciprocal overlap
                        // (i.e. sufficient overlap w.r.t. both A and B?)
                        if (!reciprocal)
                            hits.push_back(*bedItr);
                        else if (bOverlap >= overlapFraction)
                            hits.push_back(*bedItr);
                    }
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}


bool BedFile::anyHits(string chrom, CHRPOS start, CHRPOS end, string strand,
                     bool sameStrand, bool diffStrand, float overlapFraction, bool reciprocal) {

    BIN startBin, endBin;
    startBin = (start >> _binFirstShift);
    endBin = ((end-1) >> _binFirstShift);
    CHRPOS aLength = (end - start);

    /* SYNOPSIS:
    1. We loop through each UCSC BIN level for feature A's chrom.
    2. For each BIN, we loop through each B feature and return true
       if it meets all of the user's requests, which include:
       (a) overlap fractio, (b) strandedness, (c) reciprocal overlap.
       Otherwise, return false.
    */
    for (BINLEVEL i = 0; i < _binLevels; ++i) {
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {
            // move to the next bin if this one is empty
            if (bedMap[chrom][j].empty()) continue;
            vector<BED>::const_iterator bedItr = bedMap[chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[chrom][j].end();
            for (; bedItr != bedEnd; ++bedItr) {
                CHRPOS s = max(start, bedItr->start);
                CHRPOS e = min(end, bedItr->end);
                int overlapBases = (e - s); 
                // 1. is there sufficient overlap w.r.t A?
                if ( (float) overlapBases / (float) aLength  >= overlapFraction) {
                    CHRPOS bLength = (bedItr->end - bedItr->start);
                    float bOverlap = ( (float) overlapBases / (float) bLength );
                    bool strands_are_same = (strand == bedItr->strand);
                    // 2. does the overlap meet the user's strand repuirements?
                    if ( (sameStrand == false && diffStrand == false)
                        ||
                        (sameStrand == true && strands_are_same == true)
                        ||
                        (diffStrand == true && strands_are_same == false)
                        )
                    {
                        // 3. did the user request reciprocal overlap
                        // (i.e. sufficient overlap w.r.t. both A and B?)
                        if (!reciprocal)
                            return true;
                        else if (bOverlap >= overlapFraction)
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


void BedFile::countHits(const BED &a, bool sameStrand, bool diffStrand, bool countsOnly) {

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
                
                bool strands_are_same = (a.strand == bedItr->strand);
                // skip the hit if not on the same strand (and we care)
                if ((sameStrand == true && strands_are_same == false) ||
                    (diffStrand == true && strands_are_same == true)
                   ) 
                {
                    continue;
                }
                else if (overlaps(bedItr->start, bedItr->end, a.start, a.end) > 0) {

                    bedItr->count++;
                    if (countsOnly == false) {
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
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
}


void BedFile::countSplitHits(const vector<BED> &bedBlocks, bool sameStrand, bool diffStrand, bool countsOnly) {

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

                    bool strands_are_same = (blockItr->strand == bedItr->strand);
                    // skip the hit if not on the same strand (and we care)
                    if ((sameStrand == true && strands_are_same == false) ||
                        (diffStrand == true && strands_are_same == true)
                       ) 
                    {
                        continue;
                    }
                    else if (overlaps(bedItr->start, bedItr->end, blockItr->start, blockItr->end) > 0) {
                        if (countsOnly == false) {
                            if (blockItr->zeroLength == false) {
                                bedItr->depthMap[blockItr->start+1].starts++;
                                bedItr->depthMap[blockItr->end].ends++;
                            }
                            else {
                                // correct for the fact that we artificially expanded the zeroLength feature
                                bedItr->depthMap[blockItr->start+2].starts++;
                                bedItr->depthMap[blockItr->end-1].ends++;
                            }
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


void BedFile::countListHits(const BED &a, int index, bool sameStrand, bool diffStrand) {

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

                bool strands_are_same = (a.strand == bedItr->strand);
                // skip the hit if not on the same strand (and we care)
                if ((sameStrand == true && strands_are_same == false) ||
                    (diffStrand == true && strands_are_same == true)
                   ) 
                {
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


void BedFile::setBedType (int colNums) { bedType = colNums; }
void BedFile::setBed12 (bool isBed12) { this->isBed12 = isBed12; }

void BedFile::loadBedFileIntoMap() {

    BED bedEntry;
    
    Open();
    while (GetNextBed(bedEntry)) {
        if (_status == BED_VALID) {
            BIN bin = getBin(bedEntry.start, bedEntry.end);
            bedMap[bedEntry.chrom][bin].push_back(bedEntry);
        }
    }
    Close();
}


void BedFile::loadBedCovFileIntoMap() {

    BED bedEntry;
    Open();
    while (GetNextBed(bedEntry)) {
        if (_status == BED_VALID) {
            BIN bin = getBin(bedEntry.start, bedEntry.end);

            BEDCOV bedCov;
            bedCov.chrom        = bedEntry.chrom;
            bedCov.start        = bedEntry.start;
            bedCov.end          = bedEntry.end;
            bedCov.name         = bedEntry.name;
            bedCov.score        = bedEntry.score;
            bedCov.strand       = bedEntry.strand;
            bedCov.fields       = bedEntry.fields;
            bedCov.other_idxs   = bedEntry.other_idxs;
            bedCov.zeroLength   = bedEntry.zeroLength;
            bedCov.count = 0;
            bedCov.minOverlapStart = INT_MAX;

            bedCovMap[bedEntry.chrom][bin].push_back(bedCov);
        }
    }
    Close();
}

void BedFile::loadBedCovListFileIntoMap() {

    BED bedEntry;
    Open();
    while (GetNextBed(bedEntry)) {
        if (_status == BED_VALID) {
            BIN bin = getBin(bedEntry.start, bedEntry.end);

            BEDCOVLIST bedCovList;
            bedCovList.chrom        = bedEntry.chrom;
            bedCovList.start        = bedEntry.start;
            bedCovList.end          = bedEntry.end;
            bedCovList.name         = bedEntry.name;
            bedCovList.score        = bedEntry.score;
            bedCovList.strand       = bedEntry.strand;
            bedCovList.fields       = bedEntry.fields;
            bedCovList.other_idxs   = bedEntry.other_idxs;
            bedCovList.zeroLength   = bedEntry.zeroLength;

            bedCovListMap[bedEntry.chrom][bin].push_back(bedCovList);
        }
    }
    Close();
}


void BedFile::loadBedFileIntoMapNoBin() {

    BED bedEntry;
    
    Open();
    while (GetNextBed(bedEntry)) {
        if (_status == BED_VALID) {
            bedMapNoBin[bedEntry.chrom].push_back(bedEntry);
        }
    }
    Close();

    // sort the BED entries for each chromosome
    // in ascending order of start position
    for (masterBedMapNoBin::iterator m = this->bedMapNoBin.begin(); m != this->bedMapNoBin.end(); ++m) {
        sort(m->second.begin(), m->second.end(), sortByStart);
    }
}
