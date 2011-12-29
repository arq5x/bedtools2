/*****************************************************************************
  coverageMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "coverageBed.h"
#include "version.h"

using namespace std;

// define the version
#define PROGRAM_NAME "bedtools coverage"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void coverage_help(void);

int coverage_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;

    // parm flags
    bool sameStrand    = false;
    bool diffStrand    = false;
    bool writeHistogram = false;
    bool eachBase       = false;
    bool obeySplits     = false;
    bool bamInput       = false;
    bool haveBedA       = false;
    bool haveBedB       = false;
    bool countsOnly     = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) coverage_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-a", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBedA = true;
                bedAFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-abam", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveBedA = true;
                bamInput = true;
                bedAFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBedB = true;
                bedBFile = argv[i + 1];
                i++;
            }
        }
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            sameStrand = true;
        }
        else if (PARAMETER_CHECK("-S", 2, parameterLength)) {
            diffStrand = true;
        }
        else if (PARAMETER_CHECK("-hist", 5, parameterLength)) {
            writeHistogram = true;
        }
        else if(PARAMETER_CHECK("-d", 2, parameterLength)) {
            eachBase = true;
        }
        else if (PARAMETER_CHECK("-split", 6, parameterLength)) {
            obeySplits = true;
        }
        else if (PARAMETER_CHECK("-counts", 7, parameterLength)) {
            countsOnly = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBedA || !haveBedB) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -a and -b files. " << endl << "*****" << endl;
        showHelp = true;
    }
    
    if (sameStrand && diffStrand) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -s OR -S, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedCoverage *bg = new BedCoverage(bedAFile, bedBFile, sameStrand, diffStrand,
                                          writeHistogram, bamInput, obeySplits, eachBase, countsOnly);
        delete bg;
    }
    else {
        coverage_help();
    }
    return 0;
}

void coverage_help(void) {

    cerr << "\nTool:    bedtools coverage (aka coverageBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Returns the depth and breadth of coverage of features from A" << endl;
    cerr << "\t on the intervals in B." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-abam\t"         << "The A input file is in BAM format." << endl << endl;

    cerr << "\t-s\t"            << "Require same strandedness.  That is, only counts hits in A that" << endl;
    cerr                        << "\t\toverlap B on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are counted without respect to strand." << endl << endl;

    cerr << "\t-S\t"            << "Require different strandedness.  That is, only report hits in A" << endl;
    cerr                        << "\t\tthat overlap B on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are counted without respect to strand." << endl << endl;

    cerr << "\t-hist\t"         << "Report a histogram of coverage for each feature in B" << endl;
    cerr                        << "\t\tas well as a summary histogram for _all_ features in B." << endl << endl;
    cerr                        << "\t\tOutput (tab delimited) after each feature in B:" << endl;
    cerr                        << "\t\t  1) depth\n\t\t  2) # bases at depth\n\t\t  3) size of B\n\t\t  4) % of B at depth" << endl << endl;

    cerr << "\t-d\t"            << "Report the depth at each position in each B feature." << endl;
    cerr                        << "\t\tPositions reported are one based.  Each position" << endl;
    cerr                        << "\t\tand depth follow the complete B feature." << endl << endl;
    
    cerr << "\t-counts\t"       << "Only report the count of overlaps, don't compute fraction, etc." << endl << endl;

    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl;
    cerr                        << "\t\twhen computing coverage." << endl;
    cerr                        << "\t\tFor BAM files, this uses the CIGAR \"N\" and \"D\" operations " << endl;
    cerr                        << "\t\tto infer the blocks for computing coverage." << endl;
    cerr                        << "\t\tFor BED12 files, this uses the BlockCount, BlockStarts," << endl;
    cerr                        << "\t\tand BlockEnds fields (i.e., columns 10,11,12)." << endl << endl;

    cerr << "Default Output:  " << endl;
    cerr << "\t" << " After each entry in B, reports: " << endl;
    cerr << "\t   1) The number of features in A that overlapped the B interval." << endl;
    cerr << "\t   2) The number of bases in B that had non-zero coverage." << endl;
    cerr << "\t   3) The length of the entry in B." << endl;
    cerr << "\t   4) The fraction of bases in B that had non-zero coverage." << endl << endl;

    exit(1);
}
