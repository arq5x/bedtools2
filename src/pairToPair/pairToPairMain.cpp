/*****************************************************************************
  pairToPairMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "pairToPair.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools pairtopair"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void pairtopair_help(void);

int pairtopair_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;

    // input arguments
    float overlapFraction = 1E-9;
    int slop = 0;
    string searchType = "both";

    // flags to track parameters
    bool haveBedA = false;
    bool haveBedB = false;
    bool haveSearchType = false;
    bool haveFraction = false;
    bool ignoreStrand = false;
    bool requireDifferentNames = false;
    bool haveSlop = false;
    bool strandedSlop = false;
    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) pairtopair_help();

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
        else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBedB = true;
                bedBFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-type", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveSearchType = true;
                searchType = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-f", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveFraction = true;
                overlapFraction = atof(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-slop", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveSlop = true;
                slop = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-ss", 3, parameterLength)) {
            strandedSlop = true;
        }
        else if(PARAMETER_CHECK("-rdn", 4, parameterLength)) {
            requireDifferentNames = true;
        }
        else if(PARAMETER_CHECK("-is", 3, parameterLength)) {
            ignoreStrand = true;
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

    if (haveSearchType && (searchType != "neither") && (searchType != "both") && (searchType != "either") && (searchType != "notboth")) {
        cerr << endl << "*****" << endl << "*****ERROR: Request \"both\",\"neither\",\"either\",or \"notboth\"" << endl << "*****" << endl;
        showHelp = true;
    }

    if (strandedSlop == true && haveSlop == false) {
        cerr << endl << "*****" << endl << "*****ERROR: Need a -slop value if requesting -ss." << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {

        PairToPair *bi = new PairToPair(bedAFile, bedBFile, overlapFraction, searchType,
                                        ignoreStrand, requireDifferentNames, slop, strandedSlop);
        delete bi;
        return 0;
    }
    else {
        pairtopair_help();
    }
    return 0;
}


void pairtopair_help(void) {

    cerr << "\nTool:    bedtools pairtopair (aka pairToPair)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Report overlaps between two paired-end BED files (BEDPE)." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <BEDPE> -b <BEDPE>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-f\t"                    << "Minimum overlap required as fraction of A (e.g. 0.05)." << endl;
    cerr                                << "\t\tDefault is 1E-9 (effectively 1bp)." << endl << endl;

    cerr << "\t-type \t"                << "Approach to reporting overlaps between A and B." << endl << endl;
    cerr                                << "\t\tneither\tReport overlaps if neither end of A overlaps B." << endl;
    cerr                                << "\t\teither\tReport overlaps if either ends of A overlap B." << endl;
    cerr                                << "\t\tboth\tReport overlaps if both ends of A overlap B." << endl;
    cerr                                << "\t\tnotboth\tReport overlaps if one or neither of A's overlap B." << endl;
    cerr                                << "\t\t- Default = both." << endl << endl;

    cerr << "\t-slop \t"                << "The amount of slop (in b.p.). to be added to each footprint of A." << endl;
    cerr                                << "\t\t*Note*: Slop is subtracted from start1 and start2" << endl;
    cerr                                << "\t\t\tand added to end1 and end2." << endl << endl;
    cerr                                << "\t\t- Default = 0." << endl << endl;
    
    cerr << "\t-ss\t"                   << "Add slop based to each BEDPE footprint based on strand." << endl;
    cerr                                << "\t\t- If strand is \"+\", slop is only added to the end coordinates." << endl;
    cerr                                << "\t\t- If strand is \"-\", slop is only added to the start coordinates." << endl;
    cerr                                << "\t\t- By default, slop is added in both directions." << endl << endl;

    cerr << "\t-is\t"                   << "Ignore strands when searching for overlaps." << endl;
    cerr                                << "\t\t- By default, strands are enforced." << endl << endl;

    cerr << "\t-rdn\t"                  << "Require the hits to have different names (i.e. avoid self-hits)." << endl;
    cerr                                << "\t\t- By default, same names are allowed." << endl << endl;


    cerr << "Refer to the BEDTools manual for BEDPE format." << endl << endl;

    // end the program here
    exit(1);
}
