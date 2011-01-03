/*****************************************************************************
  mergeMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "mergeBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "mergeBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile  = "stdin";
    int maxDistance = 0;

    // input arguments
    bool haveBed         = true;
    bool numEntries      = false;
    bool haveMaxDistance = false;
    bool forceStrand     = false;
    bool reportNames     = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) ShowHelp();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-n", 2, parameterLength)) {
            numEntries = true;
        }
        else if(PARAMETER_CHECK("-d", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveMaxDistance = true;
                maxDistance = atoi(argv[i + 1]);
                i++;
            }
        }
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            forceStrand = true;
        }
        else if (PARAMETER_CHECK("-nms", 4, parameterLength)) {
            reportNames = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBed) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i BED file. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (reportNames && numEntries) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -n OR -nms, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedMerge *bm = new BedMerge(bedFile, numEntries, maxDistance, forceStrand, reportNames);
        delete bm;
        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Merges overlapping BED/GFF/VCF entries into a single interval." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-s\t"            << "Force strandedness.  That is, only merge features" << endl;
    cerr                        << "\t\tthat are the same strand." << endl;
    cerr                        << "\t\t- By default, merging is done without respect to strand." << endl << endl;

    cerr << "\t-n\t"            << "Report the number of BED entries that were merged." << endl;
    cerr                        << "\t\t- Note: \"1\" is reported if no merging occurred." << endl << endl;


    cerr << "\t-d\t"            << "Maximum distance between features allowed for features" << endl;
    cerr                        << "\t\tto be merged." << endl;
    cerr                        << "\t\t- Def. 0. That is, overlapping & book-ended features are merged." << endl;
    cerr                        << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-nms\t"          << "Report the names of the merged features separated by semicolons." << endl << endl;


    // end the program here
    exit(1);

}
