/*****************************************************************************
  mergeBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "mergeBed.h"



void ReportMergedNames(const map<string, bool> &names) {
    unsigned int n = 0;
    printf("\t");
    map<string, bool>::const_iterator nameItr = names.begin();
    map<string, bool>::const_iterator nameEnd = names.end();
    for (; nameItr != nameEnd; ++nameItr) {
        if (n < (names.size() - 1)) {
            cout << nameItr->first << ";";
        }
        else {
            cout << nameItr->first;
        }
        n++;
    }
}

// ===============
// = Constructor =
// ===============
BedMerge::BedMerge(string &bedFile, bool &numEntries, int &maxDistance, bool &forceStrand, bool &reportNames) {

    _bedFile = bedFile;
    _numEntries = numEntries;
    _maxDistance = maxDistance;
    _forceStrand = forceStrand;
    _reportNames = reportNames;

    _bed = new BedFile(bedFile);

    if (_forceStrand == false)
        MergeBed();
    else
        MergeBedStranded();
}


// =================
// =   Destructor  =
// =================
BedMerge::~BedMerge(void) {
}


// ===============================================
// Convenience method for reporting merged blocks
// ================================================
void BedMerge::Report(string chrom, int start, int end, const map<string, bool> &names, int mergeCount) {
    if (_bed->isZeroBased == false) {start++;}
    
    printf("%s\t%d\t%d", chrom.c_str(), start, end);
    if (_numEntries == false && _reportNames == false) {
        printf("\n");
    }
    else if (_numEntries) {
        printf("\t%d\n", mergeCount);
    }
    else if (_reportNames) {
        ReportMergedNames(names);
        printf("\n");
    }
}


// =========================================================
// Convenience method for reporting merged blocks by strand
// =========================================================
void BedMerge::ReportStranded(string chrom, int start, int end, const map<string, bool> &names, int mergeCount, string strand) {
    if (_bed->isZeroBased == false) {start++;}
    
    printf("%s\t%d\t%d", chrom.c_str(), start, end);
    if (_numEntries == false && _reportNames == false) {
        printf("\t%s\n", strand.c_str());
    }
    else if (_numEntries) {
        printf("\t%d\t%s\n", mergeCount, strand.c_str());
    }
    else if (_reportNames) {
        ReportMergedNames(names);
        printf("\t%s\n", strand.c_str());
    }
}


// =====================================================
// = Merge overlapping BED entries into a single entry =
// =====================================================
void BedMerge::MergeBed() {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bed->loadBedFileIntoMapNoBin();

    // loop through each chromosome and merge their BED entries
    masterBedMapNoBin::const_iterator m    = _bed->bedMapNoBin.begin();
    masterBedMapNoBin::const_iterator mEnd = _bed->bedMapNoBin.end();
    for (; m != mEnd; ++m) {

        // bedList is already sorted by start position.
        string chrom        = m->first;
        vector<BED> bedList = m->second;
        int mergeCount = 1;
        map<string, bool> names;

        // merge overlapping features for this chromosome.
        int start = -1;
        int end   = -1;
        vector<BED>::const_iterator bedItr = bedList.begin();
        vector<BED>::const_iterator bedEnd = bedList.end();
        for (; bedItr != bedEnd; ++bedItr) {
            // new block, no overlap
            if ( (((int) bedItr->start - end) > _maxDistance) || (end < 0)) {
                if (start >= 0) {
                    Report(chrom, start, end, names, mergeCount);
                    // reset
                    mergeCount = 1;
                    names.clear();
                }
                start = bedItr->start;
                end   = bedItr->end;
                names[bedItr->name] = true;
            }
            // same block, overlaps
            else {
                if ((int) bedItr-> end > end) end = bedItr->end;
                mergeCount++;
                names[bedItr->name] = true;
            }
        }
        if (start >= 0) {
            Report(chrom, start, end, names, mergeCount);
        }
    }
}


// ==================================================================================
// = Merge overlapping BED entries into a single entry, accounting for strandedness =
// ==================================================================================
void BedMerge::MergeBedStranded() {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bed->loadBedFileIntoMapNoBin();

    // loop through each chromosome and merge their BED entries
    masterBedMapNoBin::const_iterator m    = _bed->bedMapNoBin.begin();
    masterBedMapNoBin::const_iterator mEnd = _bed->bedMapNoBin.end();
    for (; m != mEnd; ++m) {
        
        // bedList is already sorted by start position.
        string chrom        = m->first;
        vector<BED> bedList = m->second;

        // make a list of the two strands to merge separately.
        vector<string> strands(2);
        strands[0] = "+";
        strands[1] = "-";

        // do two passes, one for each strand.
        for (unsigned int s = 0; s < strands.size(); s++) {

            int mergeCount = 1;
            int numOnStrand = 0;
            map<string, bool> names;

            // merge overlapping features for this chromosome.
            int start = -1;
            int end   = -1;
            vector<BED>::const_iterator bedItr = bedList.begin();
            vector<BED>::const_iterator bedEnd = bedList.end();
            for (; bedItr != bedEnd; ++bedItr) {

                // if forcing strandedness, move on if the hit
                // is not on the current strand.
                if (bedItr->strand != strands[s]) { continue; }
                else { numOnStrand++; }
                
                if ( (((int) bedItr->start - end) > _maxDistance) || (end < 0)) {
                    if (start >= 0) {
                        ReportStranded(chrom, start, end, names, mergeCount, strands[s]);
                        // reset
                        mergeCount = 1;
                        names.clear();
                    }
                    start = bedItr->start;
                    end   = bedItr->end;
                    names[bedItr->name] = true;
                }
                else {
                    if ((int) bedItr-> end > end) end = bedItr->end;
                    mergeCount++;
                    names[bedItr->name] = true;
                }
            }
            if (start >= 0) {
                ReportStranded(chrom, start, end, names, mergeCount, strands[s]);
            }
        }
    }
}
