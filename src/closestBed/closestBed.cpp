/*****************************************************************************
  closestBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "closestBed.h"

const int MAXSLOP = 256000000;  // 2*MAXSLOP = 512 megabases.
                                // We don't want to keep looking if we
                                // can't find a nearby feature within 512 Mb.
const int SLOPGROWTH = 2048000;


/*
    Constructor
*/
BedClosest::BedClosest(string &bedAFile, string &bedBFile, bool forceStrand, string &tieMode, bool reportDistance) {

    _bedAFile       = bedAFile;
    _bedBFile       = bedBFile;
    _forceStrand    = forceStrand;
    _tieMode        = tieMode;
    _reportDistance = reportDistance;

    _bedA           = new BedFile(bedAFile);
    _bedB           = new BedFile(bedBFile);

    FindClosestBed();
}

/*
    Destructor
*/
BedClosest::~BedClosest(void) {
}


void BedClosest::FindWindowOverlaps(BED &a, vector<BED> &hits) {

    int slop = 0;  // start out just looking for overlaps
                   // within the current bin (~128Kb)

    // update the current feature's start and end

    CHRPOS aFudgeStart = 0;
    CHRPOS aFudgeEnd;
    int numOverlaps = 0;
    vector<BED> closestB;
    float maxOverlap = 0;
    CHRPOS minDistance = INT_MAX;
    vector<CHRPOS> distances;

    // is there at least one feature in B on the same chrom
    // as the current A feature?
    if(_bedB->bedMap.find(a.chrom) != _bedB->bedMap.end()) {

        while ((numOverlaps == 0) && (slop <= MAXSLOP)) {

            // add some slop (starting at 0 bases) to a in hopes
            // of finding a hit in B
            if ((static_cast<int>(a.start) - slop) > 0)
                aFudgeStart = a.start - slop;
            else
                aFudgeStart = 0;

            if ((static_cast<int>(a.start) + slop) < (2 * MAXSLOP))
                aFudgeEnd = a.end + slop;
            else
                aFudgeEnd = 2 * MAXSLOP;

            // THE HEAVY LIFTING
            // search for hits with the current slop added
            _bedB->FindOverlapsPerBin(a.chrom, aFudgeStart, aFudgeEnd, a.strand, hits, _forceStrand);

            vector<BED>::const_iterator h = hits.begin();
            vector<BED>::const_iterator hitsEnd = hits.end();
            for (; h != hitsEnd; ++h) {

                numOverlaps++;

                // do the actual features overlap?
                int s = max(a.start, h->start);
                int e = min(a.end, h->end);
                int overlapBases = (e - s);             // the number of overlapping bases b/w a and b
                int aLength = (a.end - a.start);        // the length of a in b.p.

                // there is overlap
                if (s < e) {
                    // is there enough overlap (default ~ 1bp)
                    float overlap = (float) overlapBases / (float) aLength;
                    if ( overlap > 0 ) {
                        // is this hit the closest?
                        if (overlap > maxOverlap) {
                            maxOverlap = overlap;

                            closestB.clear();
                            closestB.push_back(*h);
                            distances.clear();
                            distances.push_back(0);
                        }
                        else if (overlap == maxOverlap) {
                            closestB.push_back(*h);
                            distances.push_back(0);
                        }
                    }
                }
                // the hit is to the "left" of A
                else if (h->end < a.start){
                    if ((a.start - h->end) < minDistance) {
                        minDistance = a.start - h->end;

                        closestB.clear();
                        closestB.push_back(*h);
                        distances.clear();
                        distances.push_back(minDistance);
                    }
                    else if ((a.start - h->end) == minDistance) {
                        closestB.push_back(*h);
                        distances.push_back(minDistance);
                    }
                }
                // the hit is to the "right" of A
                else {
                    if ((h->start - a.end) < minDistance) {
                        minDistance = h->start - a.end;

                        closestB.clear();
                        closestB.push_back(*h);
                        distances.clear();
                        distances.push_back(minDistance);
                    }
                    else if ((h->start - a.end) == minDistance) {
                        closestB.push_back(*h);
                        distances.push_back(minDistance);
                    }
                }
            }
            // if no overlaps were found, we'll widen the range
            // by SLOPGROWTH in each direction and search again.
            slop += SLOPGROWTH;
        }
    }
    // there is no feature in B on the same chromosome as A
    else {
        _bedA->reportBedTab(a);
        if (_reportDistance == true) {
            _bedB->reportNullBedTab();
            cout << -1 << endl;
        }
        else
            _bedB->reportNullBedNewLine();
    }

    // report the closest feature(s) in B to the current A feature.
    // obey the user's reporting request (_tieMode)
    if (numOverlaps > 0) {

        if (closestB.size() == 1 || _tieMode == "first") {
            _bedA->reportBedTab(a);
            if (_reportDistance == true) {
                _bedB->reportBedTab(closestB[0]);
                cout << distances[0] << endl;
            }
            else
                _bedB->reportBedNewLine(closestB[0]);
        }
        else {
            if (_tieMode == "all") {
                size_t i = 0;
                for (vector<BED>::iterator b = closestB.begin(); b != closestB.end(); ++b) {
                    _bedA->reportBedTab(a);
                    if (_reportDistance == true) {
                        _bedB->reportBedTab(*b);
                        cout << distances[i++] <<endl;
                    }
                    else
                        _bedB->reportBedNewLine(*b);
                }
            }
            else if (_tieMode == "last") {
                _bedA->reportBedTab(a);
                if (_reportDistance == true) {
                    _bedB->reportBedTab(closestB[closestB.size()-1]);
                    cout << distances[distances.size() - 1]<<endl;
                }
                else
                    _bedB->reportBedNewLine(closestB[closestB.size()-1]);
            }
        }
    }
}


void BedClosest::FindClosestBed() {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bedB->loadBedFileIntoMap();

    BED a, nullBed;
    int lineNum = 0;                    // current input line number
    vector<BED> hits;                   // vector of potential hits
    hits.reserve(100);
    BedLineStatus bedStatus;

    _bedA->Open();
    // process each entry in A in search of the closest feature in B
    while ((bedStatus = _bedA->GetNextBed(a, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            FindWindowOverlaps(a, hits);
            hits.clear();
            a = nullBed;
        }
    }
    _bedA->Close();
}
// END ClosestBed

