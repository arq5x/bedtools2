/*****************************************************************************
  closestMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "closestBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "closestBed"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;
    string tieMode = "all";

    bool haveBedA       = false;
    bool haveBedB       = false;
    bool haveTieMode    = false;
    bool forceStrand    = false;
    bool reportDistance = false;


    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if( (PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) ShowHelp();

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
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            forceStrand = true;
        }
        else if (PARAMETER_CHECK("-d", 2, parameterLength)) {
            reportDistance = true;
        }
        else if (PARAMETER_CHECK("-t", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveTieMode = true;
                tieMode = argv[i + 1];
                i++;
            }
        }
    }

    // make sure we have both input files
    if (!haveBedA || !haveBedB) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -a and -b files. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (haveTieMode && (tieMode != "all") && (tieMode != "first")
                    && (tieMode != "last")) {
        cerr << endl << "*****" << endl << "*****ERROR: Request \"all\" or \"first\" or \"last\" for Tie Mode (-t)" << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedClosest *bc = new BedClosest(bedAFile, bedBFile, forceStrand, tieMode, reportDistance);
        delete bc;
        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Authors: Aaron Quinlan (aaronquinlan@gmail.com)" << endl;
    cerr       << "\t Erik Arner, Riken" << endl << endl;

    cerr << "Summary: For each feature in A, finds the closest " << endl;
    cerr << "\t feature (upstream or downstream) in B." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-s\t"            << "Force strandedness.  That is, find the closest feature in B" << endl;
    cerr                        << "\t\tthat overlaps A on the same strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-d\t"            << "In addition to the closest feature in B, " << endl;
    cerr                        << "\t\treport its distance to A as an extra column." << endl;
    cerr                        << "\t\t- The reported distance for overlapping features will be 0." << endl << endl;


    cerr << "\t-t\t"            << "How ties for closest feature are handled.  This occurs when two" << endl;
    cerr                        << "\t\tfeatures in B have exactly the same overlap with A." << endl;
    cerr                        << "\t\tBy default, all such features in B are reported." << endl;
    cerr                        << "\t\tHere are all the options:" << endl;
    cerr                        << "\t\t- \"all\"  Report all ties (default)." << endl;
    cerr                        << "\t\t- \"first\"  Report the first tie that occurred in the B file." << endl;
    cerr                        << "\t\t- \"last\"  Report the last tie that occurred in the B file." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\tReports \"none\" for chrom and \"-1\" for all other fields when a feature" << endl;
    cerr << "\tis not found in B on the same chromosome as the feature in A." << endl;
    cerr << "\tE.g. none\t-1\t-1" << endl << endl;

    // end the program here
    exit(1);
}
