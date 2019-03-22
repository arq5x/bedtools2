/*****************************************************************************
  pairToPair.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "pairToPair.h"


/*
    Constructor
*/
PairToPair::PairToPair(string &bedAFilePE, string &bedBFilePE, float &overlapFraction,
                           string searchType, bool ignoreStrand, bool reqDiffNames, int slop, bool strandedSlop) {

    _bedAFilePE      = bedAFilePE;
    _bedBFilePE      = bedBFilePE;
    _overlapFraction = overlapFraction;
    _searchType      = searchType;
    _ignoreStrand    = ignoreStrand;
    _reqDiffNames    = reqDiffNames;
    _slop            = slop;
    _strandedSlop    = strandedSlop;

    _bedA = new BedFilePE(bedAFilePE);
    _bedB = new BedFilePE(bedBFilePE);

    IntersectPairs();
}


/*
    Destructor
*/
PairToPair::~PairToPair(void) {
}



void PairToPair::IntersectPairs() {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bedB->loadBedPEFileIntoMap();

    int lineNum = 0;
    BedLineStatus bedStatus;
    BEDPE a, nullBedPE;

    _bedA->Open();
    while ((bedStatus = _bedA->GetNextBedPE(a, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {
            // identify overlaps b/w the pairs
            FindOverlaps(a);
            a = nullBedPE;
        }
    }
    _bedA->Close();
}
// END IntersectPE



void PairToPair::FindOverlaps(const BEDPE &a) {
    //
    vector<MATE> hitsA1B1, hitsA1B2, hitsA2B1, hitsA2B2;

    // add the appropriate slop to the starts and ends
    CHRPOS start1 = a.start1;
    CHRPOS end1   = a.end1;
    CHRPOS start2 = a.start2;
    CHRPOS end2   = a.end2;

    if (_strandedSlop == true) {
        if (a.strand1 == "+")
            end1   += _slop;
        else
            start1 -= _slop;
        if (a.strand2 == "+")
            end2   += _slop;
        else
            start2 -= _slop;
    }
    else {
        (start1 - _slop) >= 0 ? start1 -= _slop : start1 = 0;
        (start2 - _slop) >= 0 ? start2 -= _slop : start2 = 0;
        end1   += _slop;
        end2   += _slop;
    }

    // Find the _potential_ hits between each end of A and B
    _bedB->FindOverlapsPerBin(1, a.chrom1, start1, end1, a.name, a.strand1, hitsA1B1, _overlapFraction, !(_ignoreStrand), _reqDiffNames);   // hits b/w A1 & B1
    _bedB->FindOverlapsPerBin(1, a.chrom2, start2, end2, a.name, a.strand2, hitsA2B1, _overlapFraction, !(_ignoreStrand), _reqDiffNames);   // hits b/w A2 & B1
    _bedB->FindOverlapsPerBin(2, a.chrom1, start1, end1, a.name, a.strand1, hitsA1B2, _overlapFraction, !(_ignoreStrand), _reqDiffNames);   // hits b/w A1 & B2
    _bedB->FindOverlapsPerBin(2, a.chrom2, start2, end2, a.name, a.strand2, hitsA2B2, _overlapFraction, !(_ignoreStrand), _reqDiffNames);   // hits b/w A2 & B2

    CHRPOS matchCount1 = (hitsA1B1.size() + hitsA2B2.size());
    CHRPOS matchCount2 = (hitsA2B1.size() + hitsA1B2.size());

    
    // report the fact that no hits were found iff _searchType is neither.
    if ((matchCount1 == 0) && (matchCount2 == 0) && (_searchType == "neither")) {
        _bedA->reportBedPENewLine(a);
    }
    else if (_searchType == "both")  {
        bool found1 = false;
        bool found2 = false;
        if ((hitsA1B1.size() > 0) || (hitsA2B2.size() > 0))
            found1 = FindHitsOnBothEnds(a, hitsA1B1, hitsA2B2);
        if ((hitsA2B1.size() > 0) || (hitsA1B2.size() > 0))
            found2 = FindHitsOnBothEnds(a, hitsA2B1, hitsA1B2);
    }
    else if (_searchType == "notboth")  {
        bool found1 = false;
        bool found2 = false;
        if ((hitsA1B1.size() > 0) || (hitsA2B2.size() > 0))
            found1 = FindHitsOnBothEnds(a, hitsA1B1, hitsA2B2);
        if ((hitsA2B1.size() > 0) || (hitsA1B2.size() > 0))
            found2 = FindHitsOnBothEnds(a, hitsA2B1, hitsA1B2);
        if (found1 == false && found2 == false)
            _bedA->reportBedPENewLine(a);
    }
    else if (_searchType == "either") {
        FindHitsOnEitherEnd(a, hitsA1B1, hitsA2B2);
        FindHitsOnEitherEnd(a, hitsA2B1, hitsA1B2);
    }
}


bool PairToPair::FindHitsOnBothEnds(const BEDPE &a, const vector<MATE> &qualityHitsEnd1,
                                    const vector<MATE> &qualityHitsEnd2) {

    map<unsigned int, vector<MATE>, less<unsigned int> > hitsMap;

    for (vector<MATE>::const_iterator h = qualityHitsEnd1.begin(); h != qualityHitsEnd1.end(); ++h) {
        hitsMap[h->lineNum].push_back(*h);
    }
    for (vector<MATE>::const_iterator h = qualityHitsEnd2.begin(); h != qualityHitsEnd2.end(); ++h) {
        hitsMap[h->lineNum].push_back(*h);
    }


    bool bothFound = false;
    for (map<unsigned int, vector<MATE>, less<unsigned int> >::iterator m = hitsMap.begin(); m != hitsMap.end(); ++m) {
        
        // hits on both sides
        if (m->second.size() >= 2) {
            bothFound = true;
            MATE b1 = m->second[0];
            MATE b2 = m->second[1];

            if (_searchType == "both") {
                _bedA->reportBedPETab(a);
                printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s",
                    b1.bed.chrom.c_str(), b1.bed.start, b1.bed.end,
                    b2.bed.chrom.c_str(), b2.bed.start, b2.bed.end,
                    b1.bed.name.c_str(), b1.bed.score.c_str(),
                    b1.bed.strand.c_str(), b2.bed.strand.c_str());
                for (size_t i = 0; i < b1.bed.other_idxs.size(); ++i)
                    printf("\t%s", b1.bed.fields[b1.bed.other_idxs[i]].c_str());
                printf("\n");
            }
        }
    }
    return bothFound;
}


void PairToPair::FindHitsOnEitherEnd(const BEDPE &a, const vector<MATE> &qualityHitsEnd1,
                                    const vector<MATE> &qualityHitsEnd2) {

    map<unsigned int, vector<MATE>, less<unsigned int> > hitsMap;

    for (vector<MATE>::const_iterator h = qualityHitsEnd1.begin(); h != qualityHitsEnd1.end(); ++h) {
        hitsMap[h->lineNum].push_back(*h);
    }
    for (vector<MATE>::const_iterator h = qualityHitsEnd2.begin(); h != qualityHitsEnd2.end(); ++h) {
        hitsMap[h->lineNum].push_back(*h);
    }

    for (map<unsigned int, vector<MATE>, less<unsigned int> >::iterator m = hitsMap.begin(); m != hitsMap.end(); ++m) {
        if (m->second.size() >= 1) {

            if ((m->second.size()) == 2) {
                MATE b1 = m->second[0];
                MATE b2 = m->second[1];

                _bedA->reportBedPETab(a);
                printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s",
                    b1.bed.chrom.c_str(), b1.bed.start, b1.bed.end,
                    b2.bed.chrom.c_str(), b2.bed.start, b2.bed.end,
                    b1.bed.name.c_str(), b1.bed.score.c_str(),
                    b1.bed.strand.c_str(), b2.bed.strand.c_str());
                for (size_t i = 0; i < b1.bed.other_idxs.size(); ++i)
                    printf("\t%s", b1.bed.fields[b1.bed.other_idxs[i]].c_str());
                printf("\n");
            }
            else {
                MATE b1 = m->second[0];

                _bedA->reportBedPETab(a);
                printf("%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%" PRId_CHRPOS "\t%" PRId_CHRPOS "\t%s\t%s\t%s\t%s",
                    b1.bed.chrom.c_str(), b1.bed.start, b1.bed.end,
                    b1.mate->bed.chrom.c_str(), b1.mate->bed.start, b1.mate->bed.end,
                    b1.bed.name.c_str(), b1.bed.score.c_str(),
                    b1.bed.strand.c_str(), b1.mate->bed.strand.c_str());
                for (size_t i = 0; i < b1.bed.other_idxs.size(); ++i)
                    printf("\t%s", b1.bed.fields[b1.bed.other_idxs[i]].c_str());
                printf("\n");
            }
        }
    }
}
