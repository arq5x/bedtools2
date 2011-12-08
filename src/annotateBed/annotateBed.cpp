/*****************************************************************************
  annotateBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "annotateBed.h"

// build
BedAnnotate::BedAnnotate(const string &mainFile, const vector<string> &annoFileNames,
            const vector<string> &annoTitles, bool sameStrand, bool diffStrand, bool reportCounts, bool reportBoth) :

    _mainFile(mainFile),
    _annoFileNames(annoFileNames),
    _annoTitles(annoTitles),
    _sameStrand(sameStrand),
    _diffStrand(diffStrand),
    _reportCounts(reportCounts),
    _reportBoth(reportBoth)
{
    _bed = new BedFile(_mainFile);
}


// destroy and delete the open file pointers
BedAnnotate::~BedAnnotate(void) {
    delete _bed;
    CloseAnnoFiles();
}


void BedAnnotate::OpenAnnoFiles() {
    for (size_t i=0; i < _annoFileNames.size(); ++i) {
        BedFile *file = new BedFile(_annoFileNames[i]);
        file->Open();
        _annoFiles.push_back(file);
    }
}


void BedAnnotate::CloseAnnoFiles() {
    for (size_t i=0; i < _annoFiles.size(); ++i) {
        BedFile *file = _annoFiles[i];
        delete file;
        _annoFiles[i] = NULL;
    }
}


void BedAnnotate::PrintHeader() {
    // print a hash to indicate header and then write a tab
    // for each field in the main file.
    printf("#");
    for (size_t i = 0; i < _bed->bedType; ++i)
        printf("\t");

    // now print the label for each file.
    if (_reportBoth == false) {
        for (size_t i = 0; i < _annoTitles.size(); ++i)
            printf("%s\t", _annoTitles[i].c_str());
        printf("\n");
    }
    else {
        for (size_t i = 0; i < _annoTitles.size(); ++i)
            printf("%s_cnt\t%s_pct", _annoTitles[i].c_str(), _annoTitles[i].c_str());
        printf("\n");
    }
}


void BedAnnotate::InitializeMainFile() {
    // process each chromosome
    masterBedCovListMap::iterator chromItr = _bed->bedCovListMap.begin();
    masterBedCovListMap::iterator chromEnd = _bed->bedCovListMap.end();
    for (; chromItr != chromEnd; ++chromItr) {
        // for each chrom, process each bin
        binsToBedCovLists::iterator binItr = chromItr->second.begin();
        binsToBedCovLists::iterator binEnd = chromItr->second.end();
        for (; binItr != binEnd; ++binItr) {
            // initialize BEDCOVLIST in this chrom/bin
            vector<BEDCOVLIST>::iterator bedItr = binItr->second.begin();
            vector<BEDCOVLIST>::iterator bedEnd = binItr->second.end();
            for (; bedItr != bedEnd; ++bedItr) {
                // initialize the depthMaps, counts, etc. for each anno file.
                for (size_t i = 0; i < _annoFiles.size(); ++i) {
                    map<unsigned int, DEPTH> dummy;
                    bedItr->depthMapList.push_back(dummy);
                    bedItr->counts.push_back(0);
                    bedItr->minOverlapStarts.push_back(INT_MAX);
                }
            }
        }
    }
}


void BedAnnotate::AnnotateBed() {

    // load the "main" bed file into a map so
    // that we can easily compare each annoFile to it for overlaps
    _bed->loadBedCovListFileIntoMap();
    // open the annotations files for processing;
    OpenAnnoFiles();
    // initialize counters, depths, etc. for the main file
    InitializeMainFile();

    // annotate the main file with the coverage from the annotation files.
    for (size_t annoIndex = 0; annoIndex < _annoFiles.size(); ++annoIndex) {
        // grab the current annotation file.
        BedFile *anno = _annoFiles[annoIndex];
        BED a;
        // process each entry in the current anno file
        while (anno->GetNextBed(a)) {
            if (anno->_status == BED_VALID) {
                _bed->countListHits(a, annoIndex, _sameStrand, _diffStrand);
            }
        }
    }

    // report the annotations of the main file from the anno file.
    ReportAnnotations();
    // close the annotations files;
    CloseAnnoFiles();
}


void BedAnnotate::ReportAnnotations() {

    if (_annoTitles.size() > 0) {
        PrintHeader();
    }

    // process each chromosome
    masterBedCovListMap::const_iterator chromItr = _bed->bedCovListMap.begin();
    masterBedCovListMap::const_iterator chromEnd = _bed->bedCovListMap.end();
    for (; chromItr != chromEnd; ++chromItr) {
        // for each chrom, process each bin
        binsToBedCovLists::const_iterator binItr = chromItr->second.begin();
        binsToBedCovLists::const_iterator binEnd = chromItr->second.end();
        for (; binItr != binEnd; ++binItr) {
            // for each chrom & bin, compute and report
            // the observed coverage for each feature
            vector<BEDCOVLIST>::const_iterator bedItr = binItr->second.begin();
            vector<BEDCOVLIST>::const_iterator bedEnd = binItr->second.end();
            for (; bedItr != bedEnd; ++bedItr) {
                // print the main BED entry.
                _bed->reportBedTab(*bedItr);

                // now report the coverage from each annotation file.
                for (size_t i = 0; i < _annoFiles.size(); ++i) {
                    unsigned int totalLength = 0;
                    int zeroDepthCount = 0; // number of bases with zero depth
                    int depth          = 0; // tracks the depth at the current base

                    // the start is either the first base in the feature OR
                    // the leftmost position of an overlapping feature. e.g. (s = start):
                    // A    ----------
                    // B    s    ------------
                    int start          = min(bedItr->minOverlapStarts[i], bedItr->start);

                    map<unsigned int, DEPTH>::const_iterator depthItr;
                    map<unsigned int, DEPTH>::const_iterator depthEnd;

                    // compute the coverage observed at each base in the feature marching from start to end.
                    for (CHRPOS pos = start+1; pos <= bedItr->end; pos++) {
                        // map pointer grabbing the starts and ends observed at this position
                        depthItr = bedItr->depthMapList[i].find(pos);
                        depthEnd = bedItr->depthMapList[i].end();

                        // increment coverage if starts observed at this position.
                        if (depthItr != depthEnd)
                            depth += depthItr->second.starts;
                        // update zero depth
                        if ((pos > bedItr->start) && (pos <= bedItr->end) && (depth == 0))
                            zeroDepthCount++;
                        // decrement coverage if ends observed at this position.
                        if (depthItr != depthEnd)
                            depth = depth - depthItr->second.ends;
                    }
                    // Summarize the coverage for the current interval,
                    CHRPOS length     = bedItr->end - bedItr->start;
                    totalLength       += length;
                    int nonZeroBases   = (length - zeroDepthCount);
                    float fractCovered = (float) nonZeroBases / length;
                    if (_reportCounts == false && _reportBoth == false)
                        printf("%f\t", fractCovered);
                    else if (_reportCounts == true && _reportBoth == false)
                        printf("%d\t", bedItr->counts[i]);
                    else if (_reportCounts == false && _reportBoth == true)
                        printf("%d\t%f\t", bedItr->counts[i], fractCovered);
                }
                // print newline for next feature.
                printf("\n");
            }
        }
    }
}


