/*****************************************************************************
  subtractMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "subtractBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools subtract"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void subtract_help(void);

int subtract_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;

    // input arguments
    float overlapFraction = 1E-9;

    bool haveBedA = false;
    bool haveBedB = false;
    bool haveFraction = false;
    bool sameStrand = false;
    bool diffStrand = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) subtract_help();

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
        else if(PARAMETER_CHECK("-f", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveFraction = true;
                overlapFraction = atof(argv[i + 1]);
                i++;
            }
        }
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            sameStrand = true;
        }
        else if (PARAMETER_CHECK("-S", 2, parameterLength)) {
            diffStrand = true;
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

        BedSubtract *bs = new BedSubtract(bedAFile, bedBFile, overlapFraction, sameStrand, diffStrand);
        delete bs;
        return 0;
    }
    else {
        subtract_help();
    }
    return 0;
}

void subtract_help(void) {

    cerr << "\nTool:    bedtools subtract (aka subtractBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Removes the portion(s) of an interval that is overlapped" << endl;
    cerr << "\t by another feature(s)." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of A." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- (FLOAT) (e.g. 0.50)" << endl << endl;

    cerr << "\t-s\t"            << "Require same strandedness.  That is, only subtract hits in B" << endl;
    cerr                        << "\t\tthat overlap A on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are subtracted without respect to strand." << endl << endl;

    cerr << "\t-S\t"            << "Force strandedness.  That is, only subtract hits in B that" << endl;
    cerr                        << "\t\toverlap A on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are subtracted without respect to strand." << endl << endl;

    // end the program here
    exit(1);
}
