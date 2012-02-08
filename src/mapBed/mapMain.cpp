/*****************************************************************************
  mapMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "mapBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools map"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void map_help(void);

int map_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;
    int column = 5;
    string operation = "sum";
    string nullValue = ".";

    // input arguments
    float overlapFraction = 1E-9;

    bool haveBedA           = false;
    bool haveBedB           = false;
    bool haveColumn         = false;
    bool haveOperation      = false;
    bool haveFraction       = false;
    bool reciprocalFraction = false;
    bool sameStrand         = false;
    bool diffStrand         = false;
    bool printHeader        = false;
    bool choseNullValue     = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) map_help();

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
        else if(PARAMETER_CHECK("-c", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveColumn = true;
                column = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-o", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveOperation = true;
                operation = argv[i + 1];
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
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            sameStrand = true;
        }
        else if (PARAMETER_CHECK("-S", 2, parameterLength)) {
            diffStrand = true;
        }
        else if (PARAMETER_CHECK("-null", 5, parameterLength)) {
            nullValue = argv[i + 1];
            choseNullValue = true;
            i++;
        }
        else if(PARAMETER_CHECK("-header", 7, parameterLength)) {
            printHeader = true;
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

    if (reciprocalFraction && !haveFraction) {
        cerr << endl << "*****" << endl << "*****ERROR: If using -r, you need to define -f." << endl << "*****" << endl;
        showHelp = true;
    }

    if (sameStrand && diffStrand) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -s OR -S, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {

        BedMap *bm = new BedMap(bedAFile, bedBFile, column, operation,
                                       overlapFraction, sameStrand,
                                       diffStrand, reciprocalFraction,
                                       choseNullValue, nullValue, 
                                       printHeader);
        delete bm;
        return 0;
    }
    else {
        map_help();
        return 0;
    }
}

void map_help(void) {

    cerr << "\nTool:    bedtools map (aka mapBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Apply a function to a column from B intervals that overlap A." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-c\t"             << "Specify the column from the B file to map onto intervals in A." << endl;
    cerr                         << "\t\t - Default = 5." << endl << endl;

    cerr << "\t-o\t"             << "Specify the operation that should be applied to -c." << endl;
    cerr                         << "\t\t Valid operations:" << endl;
    cerr                         << "\t\t    sum, min, max," << endl;
    cerr                         << "\t\t    mean, median," << endl;
    cerr                         << "\t\t    collapse (i.e., print a comma separated list (duplicates allowed)), " << endl;
    cerr                         << "\t\t    distinct (i.e., print a comma separated list (NO duplicates allowed)), " << endl;
    cerr                         << "\t\t    count" << endl;
    cerr                         << "\t\t    count_distinct (i.e., a count of the unique values in the column), " << endl;
    cerr                         << "\t\t- Default: sum" << endl << endl;

    cerr << "\t-f\t"             << "Minimum overlap required as a fraction of A." << endl;
    cerr                         << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                         << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;
                                 
    cerr << "\t-r\t"             << "Require that the fraction overlap be reciprocal for A and B." << endl;
    cerr                         << "\t\t- In other words, if -f is 0.90 and -r is used, this requires" << endl;
    cerr                         << "\t\t  that B overlap 90% of A and A _also_ overlaps 90% of B." << endl << endl;
                                 
    cerr << "\t-s\t"             << "Require same strandedness.  That is, only report hits in B" << endl;
    cerr                         << "\t\tthat overlap A on the _same_ strand." << endl;
    cerr                         << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;
                                 
    cerr << "\t-S\t"             << "Require different strandedness.  That is, only report hits in B" << endl;
    cerr                         << "\t\tthat overlap A on the _opposite_ strand." << endl;
    cerr                         << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-null\t"          << "The value to print if no overlaps are found for an A interval." << endl;
    cerr                         << "\t\t- Default - \".\"" << endl << endl;

    cerr << "\t-header\t"        << "Print the header from the A file prior to results." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1) Both input files must be sorted by chrom, then start." << endl << endl;
    
    // end the program here
    exit(1);

}
