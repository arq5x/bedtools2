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



void BedMerge::ReportMergedNames(const vector<string> &names) {
    if (names.size() > 0) {
        printf("\t");
        vector<string>::const_iterator nameItr = names.begin();
        vector<string>::const_iterator nameEnd = names.end();
        for (; nameItr != nameEnd; ++nameItr) {
            if (nameItr < (nameEnd - 1))
                cout << *nameItr << ";";
            else
                cout << *nameItr;
        }
    }
    else {
        cerr << endl 
             << "*****" << endl 
             << "*****ERROR: No names found to report for the -names option. Exiting." << endl 
             << "*****" << endl;
        exit(1);
    }
}


void BedMerge::ReportMergedScores(const vector<string> &scores) {
    
    // setup a VectorOps instances for the list of scores.
    // VectorOps methods used for each possible operation.
    VectorOps vo(scores);
    std::stringstream buffer;
    if (scores.size() > 0) {
        if (_scoreOp == "sum")
            buffer << setprecision (PRECISION) << vo.GetSum();
        else if (_scoreOp == "min")
            buffer << setprecision (PRECISION) << vo.GetMin();
        else if (_scoreOp == "max")
            buffer << setprecision (PRECISION) << vo.GetMax();
        else if (_scoreOp == "mean")
            buffer << setprecision (PRECISION) << vo.GetMean();
        else if (_scoreOp == "median")
            buffer << setprecision (PRECISION) << vo.GetMedian();
        else if (_scoreOp == "mode")
            buffer << setprecision (PRECISION) << vo.GetMode();
        else if (_scoreOp == "antimode")
            buffer << setprecision (PRECISION) << vo.GetAntiMode();
        else if (_scoreOp == "collapse")
            buffer << setprecision (PRECISION) << vo.GetCollapse();
        cout << "\t" << buffer.str();
    }
    else {        
        cerr << endl 
             << "*****" << endl 
             << "*****ERROR: No scores found to report for the -scores option. Exiting." << endl 
             << "*****" << endl;
        exit(1);
    }
}

// ===============
// = Constructor =
// ===============
BedMerge::BedMerge(string &bedFile, 
                   bool numEntries, 
                   int  maxDistance, 
                   bool forceStrand, 
                   bool reportNames, 
                   bool reportScores,
                   const string &scoreOp) :
    _bedFile(bedFile),
    _numEntries(numEntries),
    _forceStrand(forceStrand),
    _reportNames(reportNames),
    _reportScores(reportScores),
    _scoreOp(scoreOp),
    _maxDistance(maxDistance)
{
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
void BedMerge::Report(string chrom, int start, int end, 
                      const vector<string> &names, const vector<string> &scores, int mergeCount) 
{
    // ARQ: removed to force all output to be zero-based, BED format, reagrdless of input type
    //if (_bed->isZeroBased == false) {start++;}
    
    printf("%s\t%d\t%d", chrom.c_str(), start, end);
    // just the merged intervals
    if (_numEntries == false && _reportNames == false && _reportScores == false) {
        printf("\n");
    }
    // merged intervals and counts    
    else if (_numEntries == true && _reportNames == false && _reportScores == false) {
        printf("\t%d\n", mergeCount);
    }
    // merged intervals, counts, and scores
    else if (_numEntries == true && _reportNames == false && _reportScores == true) {
        printf("\t%d", mergeCount);
        ReportMergedScores(scores);
        printf("\n");
    }
    // merged intervals, counts, and names
    else if (_numEntries == true && _reportNames == true && _reportScores == false) {
        ReportMergedNames(names);
        printf("\t%d\n", mergeCount);
    }
    // merged intervals, counts, names, and scores
    else if (_numEntries == true && _reportNames == true && _reportScores == true) {
        ReportMergedNames(names);
        ReportMergedScores(scores);
        printf("\t%d\n", mergeCount);
    }
    // merged intervals and names        
    else if (_numEntries == false && _reportNames == true && _reportScores == false) {
        ReportMergedNames(names);
        printf("\n");
    }
    // merged intervals and scores        
    else if (_numEntries == false && _reportNames == false && _reportScores == true) {
        ReportMergedScores(scores);
        printf("\n");
    }
    // merged intervals, names, and scores        
    else if (_numEntries == false && _reportNames == true && _reportScores == true) {
        ReportMergedNames(names);
        ReportMergedScores(scores);
        printf("\n");
    }
}


// =========================================================
// Convenience method for reporting merged blocks by strand
// =========================================================
void BedMerge::ReportStranded(string chrom, int start, int end, 
                              const vector<string> &names, const vector<string> &scores,
                              int mergeCount, string strand) 
{
    // ARQ: removed to force all output to be zero-based, BED format, reagrdless of input type
    //if (_bed->isZeroBased == false) {start++;}
    
    printf("%s\t%d\t%d", chrom.c_str(), start, end);
    // just the merged intervals
    if (_numEntries == false && _reportNames == false && _reportScores == false) {
        printf("\t%s\n", strand.c_str());
    }
    // merged intervals and counts    
    else if (_numEntries == true && _reportNames == false && _reportScores == false) {
        printf("\t%d\t%s\n", mergeCount, strand.c_str());
    }
    // merged intervals, counts, and scores
    else if (_numEntries == true && _reportNames == false && _reportScores == true) {
        printf("\t%d", mergeCount);
        ReportMergedScores(scores);
        printf("\t%s\n", strand.c_str());
    }
    // merged intervals, counts, and names
    else if (_numEntries == true && _reportNames == true && _reportScores == false) {
        ReportMergedNames(names);
        printf("\t%d\t%s", mergeCount, strand.c_str());
        printf("\n");
    }
    // merged intervals, counts, names, and scores
    else if (_numEntries == true && _reportNames == true && _reportScores == true) {
        ReportMergedNames(names);
        ReportMergedScores(scores);
        printf("\t%s\t%d", strand.c_str(), mergeCount);
        printf("\n");
    }
    // merged intervals and names        
    else if (_numEntries == false && _reportNames == true && _reportScores == false) {
        ReportMergedNames(names);
        printf("\t%s\n", strand.c_str());
    }
    // merged intervals and scores        
    else if (_numEntries == false && _reportNames == false && _reportScores == true) {
        ReportMergedScores(scores);
        printf("\t%s\n", strand.c_str());
    }
    // merged intervals, names, and scores        
    else if (_numEntries == false && _reportNames == true && _reportScores == true) {
        ReportMergedNames(names);
        ReportMergedScores(scores);
        printf("\t%s\n", strand.c_str());
    }
}


// =====================================================
// = Merge overlapping BED entries into a single entry =
// =====================================================
void BedMerge::MergeBed() {
    int mergeCount = 1;
    vector<string> names;
    vector<string> scores;
    int start = -1;
    int end   = -1;
    BED prev, curr;
    
    _bed->Open();
    while (_bed->GetNextBed(curr, true)) { // true = force sorted intervals
        if (_bed->_status != BED_VALID)
            continue;            
        // new block, no overlap
        if ( (((int) curr.start - end) > _maxDistance) || (curr.chrom != prev.chrom)) {
            if (start >= 0) {
                Report(prev.chrom, start, end, names, scores, mergeCount);
                // reset
                mergeCount = 1;
                names.clear();
                scores.clear();
            }
            start = curr.start;
            end   = curr.end;
            if (!curr.name.empty())
                names.push_back(curr.name);
            if (!curr.score.empty())
            scores.push_back(curr.score);
        }
        // same block, overlaps
        else {
            if ((int) curr.end > end) 
                end = curr.end;
            if (!curr.name.empty())
                names.push_back(curr.name);
            if (!curr.score.empty())
                scores.push_back(curr.score);
            mergeCount++;
        }
        prev = curr;
    }
    if (start >= 0) {
        Report(prev.chrom, start, end, names, scores, mergeCount);
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

        // make a list of the two strands to merge separately.
        vector<string> strands(2);
        strands[0] = "+";
        strands[1] = "-";
        // do two passes, one for each strand.
        for (unsigned int s = 0; s < strands.size(); s++) {
            int mergeCount = 1;
            int numOnStrand = 0;
            vector<string> names;
            vector<string> scores;

            // merge overlapping features for this chromosome.
            int start = -1;
            int end   = -1;
            vector<BED>::const_iterator bedItr = m->second.begin();
            vector<BED>::const_iterator bedEnd = m->second.end();
            for (; bedItr != bedEnd; ++bedItr) {
                // if forcing strandedness, move on if the hit
                // is not on the current strand.
                if (bedItr->strand != strands[s]) { continue; }
                else { numOnStrand++; }
                if ( (((int) bedItr->start - end) > _maxDistance) || (end < 0)) {
                    if (start >= 0) {
                        ReportStranded(chrom, start, end, names, scores, mergeCount, strands[s]);
                        // reset
                        mergeCount = 1;
                        names.clear();
                        scores.clear();
                    }
                    start = bedItr->start;
                    end   = bedItr->end;
                    if (!bedItr->name.empty())  names.push_back(bedItr->name);
                    if (!bedItr->score.empty()) scores.push_back(bedItr->score);
                }
                else {
                    if ((int) bedItr-> end > end) end = bedItr->end;
                    mergeCount++;
                    if (!bedItr->name.empty())  names.push_back(bedItr->name);
                    if (!bedItr->score.empty()) scores.push_back(bedItr->score);
                }
            }
            if (start >= 0) {
                ReportStranded(chrom, start, end, names, scores, mergeCount, strands[s]);
            }
        }
    }
}
