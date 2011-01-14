/*****************************************************************************
  intersectBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "fjoin.h"
#include <queue>
#include <set>

bool leftOf(const BED &a, const BED &b);


bool BedIntersect::processHits(BED &a, vector<BED> &hits) {
    // how many overlaps are there b/w the bed and the set of hits?
    int s, e, overlapBases;
    int  numOverlaps = 0;
    bool hitsFound   = false;
    int aLength      = (a.end - a.start);   // the length of a in b.p.

    // loop through the hits and report those that meet the user's criteria
    vector<BED>::const_iterator h       = hits.begin();
    vector<BED>::const_iterator hitsEnd = hits.end();
    for (; h != hitsEnd; ++h) {
        s            = max(a.start, h->start);
        e            = min(a.end, h->end);
        overlapBases = (e - s);             // the number of overlapping bases b/w a and b

        // is there enough overlap relative to the user's request? (default ~ 1bp)
        if ( ( (float) overlapBases / (float) aLength ) >= _overlapFraction ) {
            // Report the hit if the user doesn't care about reciprocal overlap between A and B.
            if (_reciprocal == false) {
                hitsFound = true;
                numOverlaps++;
                if (_printable == true)
                    ReportOverlapDetail(overlapBases, a, *h, s, e);
            }
            // we require there to be sufficient __reciprocal__ overlap
            else {
                int bLength    = (h->end - h->start);
                float bOverlap = ( (float) overlapBases / (float) bLength );
                if (bOverlap >= _overlapFraction) {
                    hitsFound = true;
                    numOverlaps++;
                    if (_printable == true)
                        ReportOverlapDetail(overlapBases, a, *h, s, e);
                }
            }
        }
    }
    // report the summary of the overlaps if requested.
    ReportOverlapSummary(a, numOverlaps);
    // were hits found for this BED feature?
    return hitsFound;
}

/*
    Constructor
*/
BedIntersect::BedIntersect(string bedAFile, string bedBFile, bool anyHit,
                           bool writeA, bool writeB, bool writeOverlap, bool writeAllOverlap,
                           float overlapFraction, bool noHit, bool writeCount, bool forceStrand,
                           bool reciprocal, bool obeySplits, bool bamInput, bool bamOutput) {

    _bedAFile            = bedAFile;
    _bedBFile            = bedBFile;
    _anyHit              = anyHit;
    _noHit               = noHit;
    _writeA              = writeA;
    _writeB              = writeB;
    _writeOverlap        = writeOverlap;
    _writeAllOverlap     = writeAllOverlap;
    _writeCount          = writeCount;
    _overlapFraction     = overlapFraction;
    _forceStrand         = forceStrand;
    _reciprocal          = reciprocal;
    _obeySplits          = obeySplits;
    _bamInput            = bamInput;
    _bamOutput           = bamOutput;

    if (_anyHit || _noHit || _writeCount)
        _printable = false;
    else
        _printable = true;

    // create new BED file objects for A and B
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);

    IntersectBed();
}


/*
    Destructor
*/
BedIntersect::~BedIntersect(void) {
}


bool leftOf(const BED &a, const BED &b) {
    return (a.end <= b.start);
}


void BedIntersect::ReportOverlapDetail(const int &overlapBases, const BED &a, const BED &b,
                                       const CHRPOS &s, const CHRPOS &e) {
    // default. simple intersection only
    if (_writeA == false && _writeB == false && _writeOverlap == false) {
        _bedA->reportBedRangeNewLine(a,s,e);
    }
    //  -wa -wbwrite the original A and B
    else if (_writeA == true && _writeB == true) {
        _bedA->reportBedTab(a);
        _bedB->reportBedNewLine(b);
    }
    // -wa write just the original A
    else if (_writeA == true) {
        _bedA->reportBedNewLine(a);
    }
    // -wb write the intersected portion of A and the original B
    else if (_writeB == true) {
        _bedA->reportBedRangeTab(a,s,e);
        _bedB->reportBedNewLine(b);
    }
    // -wo write the original A and B plus the no. of overlapping bases.
    else if (_writeOverlap == true) {
        _bedA->reportBedTab(a);
        _bedB->reportBedTab(b);
        printf("%d\n", overlapBases);
    }
}


void BedIntersect::ReportOverlapSummary(const BED &a, const int &numOverlapsFound) {
    // -u  just report the fact that there was >= 1 overlaps
    if (_anyHit && (numOverlapsFound >= 1)) {
        _bedA->reportBedNewLine(a);
    }
    // -c  report the total number of features overlapped in B
    else if (_writeCount) {
        _bedA->reportBedTab(a);
        printf("%d\n", numOverlapsFound);
    }
    // -v  report iff there were no overlaps
    else if (_noHit && (numOverlapsFound == 0)) {
        _bedA->reportBedNewLine(a);
    }
    // -wao the user wants to force the reporting of 0 overlap
    else if (_writeAllOverlap && (numOverlapsFound == 0)) {
        _bedA->reportBedTab(a);
        _bedB->reportNullBedTab();
        printf("0\n");
    }
}



void BedIntersect::Scan(BED *x, vector<BED *> *windowX, BedLineStatus xStatus,
                  const BED &y, vector<BED *> *windowY, BedLineStatus yStatus) {

    if (xStatus != BED_VALID) {
        return;
    }

    std::vector<BED *>::iterator wYIter = windowY->begin();
    while (wYIter != windowY->end()) {
        if (leftOf(*(*wYIter), *x) == true) {
            (*wYIter)->finished = true;
            wYIter = windowY->erase(wYIter);  // erase auto-increments to the next position
        }
        else if (overlaps((*wYIter)->start, (*wYIter)->end, x->start, x->end) > 0) {
            if (_lastPick == 0) {
                AddHits(x, *(*wYIter));
            }
            else {
                AddHits(*wYIter, *x);
            }
            ++wYIter;  // force incrementing
        }
    }
    if (leftOf(*x,y) == false)
        windowX->push_back(x);
    else {
        x->finished = true;
    }
    // dump the buffered results (if any)
    FlushOutputBuffer();
}


void BedIntersect::AddHits(BED *x, const BED &y) {
    if (x->overlaps.empty() == true)
        _outputBuffer.push(x);
    x->overlaps.push_back(y);
}


void BedIntersect::FlushOutputBuffer(bool final) {
    while (_outputBuffer.empty() == false)
    {
        if (final == false && _outputBuffer.front()->finished == false)
            break;

        processHits(*_outputBuffer.front(), _outputBuffer.front()->overlaps);
        // remove the finished BED entry from the heap
        delete _outputBuffer.front();
        _outputBuffer.pop();
    }
}


vector<BED*>* BedIntersect::GetWindow(const string &chrom, bool isA) {

    // iterator to test if a window for a given chrom exists.
    map<string, vector<BED*> >::iterator it;

    // grab the current window for A or B, depending on
    // the request.  if a window hasn't yet been created
    // for the requested chrom, create one.

    if (isA) {
        it = _windowA.find(chrom);
        if (it != _windowA.end()) {
            return & _windowA[chrom];
        }
        else {
            _windowA.insert(pair<string, vector<BED *> >(chrom, vector<BED *>()));
            return & _windowA[chrom];
        }
    }
    else {
        it = _windowB.find(chrom);
        if (it != _windowB.end()) {
            return & _windowB[chrom];
        }
        else {
            _windowB.insert(pair<string, vector<BED *> >(chrom, vector<BED *>()));
            return & _windowB[chrom];
        }
    }
}


void BedIntersect::ChromSwitch(const string &chrom) {

    vector<BED*>::iterator windowAIter = _windowA[chrom].begin();
    vector<BED*>::iterator windowAEnd  = _windowA[chrom].end();
    for (; windowAIter != windowAEnd; ++windowAIter)
        (*windowAIter)->finished = true;

    vector<BED*>::iterator windowBIter = _windowB[chrom].begin();
    vector<BED*>::iterator windowBEnd  = _windowB[chrom].end();
    for (; windowBIter != windowBEnd; ++windowBIter)
        (*windowBIter)->finished = true;

    FlushOutputBuffer();
}


void BedIntersect::IntersectBed() {

    int aLineNum = 0;
    int bLineNum = 0;

    // current feature from each file
    BED *a, *b, *prevA, *prevB;

    // status of the current lines
    BedLineStatus aStatus, bStatus;

    // open the files; get the first line from each
    _bedA->Open();
    _bedB->Open();

    prevA = NULL;
    prevB = NULL;
    a = new BED();
    b = new BED();
    aStatus = _bedA->GetNextBed(*a, aLineNum);
    bStatus = _bedB->GetNextBed(*b, bLineNum);

    while (aStatus != BED_INVALID || bStatus != BED_INVALID) {

        if ((a->start <= b->start) && (a->chrom == b->chrom)) {
            prevA = a;
            _lastPick = 0;
            Scan(a, GetWindow(a->chrom, true),  aStatus,
                *b, GetWindow(a->chrom, false), bStatus);

            a = new BED();
            aStatus = _bedA->GetNextBed(*a, aLineNum);
        }
        else if ((a->start > b->start) && (a->chrom == b->chrom)) {
            prevB = b;
            _lastPick = 1;
            Scan(b, GetWindow(b->chrom, false), bStatus,
                *a, GetWindow(b->chrom, true),  aStatus);

            b = new BED();
            bStatus = _bedB->GetNextBed(*b, bLineNum);
        }
        else if (a->chrom != b->chrom) {
            // A was most recently read
            if (_lastPick == 0) {
                prevB = b;
                while (b->chrom == prevA->chrom){
                    _windowB[prevA->chrom].push_back(b);
                    b = new BED();
                    bStatus = _bedB->GetNextBed(*b, bLineNum);
                }
                Scan(prevA, GetWindow(prevA->chrom, true),  aStatus,
                    *prevB, GetWindow(prevA->chrom, false),  bStatus);
            }
            // B was most recently read
            else {
                prevA = a;
                while (a->chrom == prevB->chrom) {
                    _windowA[prevB->chrom].push_back(a);
                    a = new BED();
                    aStatus = _bedA->GetNextBed(*a, aLineNum);
                }
                Scan(prevB, GetWindow(prevB->chrom, false), bStatus,
                    *prevA, GetWindow(prevB->chrom, true),  aStatus);
            }
            FlushOutputBuffer(true);
        }
        if (prevA!=NULL&&prevB!=NULL)
            //cout << prevA->chrom << " " << a->chrom << " " << a->start << " "
            //     << prevB->chrom << " " << b->chrom << " " << b->start << "\n";
        if (aStatus == BED_INVALID) a->start = INT_MAX;
        if (bStatus == BED_INVALID) b->start = INT_MAX;
    }

    // clear out the final bit of staged output
    FlushOutputBuffer(true);

    // close the files
    _bedA->Close();
    _bedB->Close();
}


