/*****************************************************************************
  jaccardMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "jaccard.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools jaccard"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void jaccard_help(void);

int jaccard_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;

    // input arguments
    float overlapFraction = 1E-9;

    bool haveBedA           = false;
    bool haveBedB           = false;
    bool haveFraction       = false;
    bool reciprocalFraction = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) jaccard_help();

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
        else if(PARAMETER_CHECK("-r", 2, parameterLength)) {
            reciprocalFraction = true;
        }
        else {
            cerr << endl 
                 << "*****ERROR: Unrecognized parameter: " 
                 << argv[i] 
                 << " *****" 
                 << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBedA || !haveBedB) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -a and -b files. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (reciprocalFraction && !haveFraction) {
        cerr << endl << "*****" << endl << "*****ERROR: If using -r, you need to define -f." << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {

        Jaccard *j = new Jaccard(bedAFile, bedBFile,
                                 overlapFraction, reciprocalFraction);
        delete j;
        return 0;
    }
    else {
        jaccard_help();
        return 0;
    }
}

void jaccard_help(void) {

    cerr << "\nTool:    bedtools jaccard (aka jaccard)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Calculate Jaccard statistic b/w two feature files."
         << endl
         << "\t Jaccard is the length of the intersection over the union."
         << endl
         << "\t Values range from 0 (no intersection) to 1 (self intersection)."
         << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;


    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of A." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-r\t"            << "Require that the fraction overlap be reciprocal for A and B." << endl;
    cerr                        << "\t\t- In other words, if -f is 0.90 and -r is used, this requires" << endl;
    cerr                        << "\t\t  that B overlap 90% of A and A _also_ overlaps 90% of B." << endl << endl;
 
    cerr << "Notes: " << endl;
    cerr << "\t(1) Input files must be sorted by chrom, then start position." 
         << endl << endl;

    // end the program here
    exit(1);

}
