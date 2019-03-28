/*****************************************************************************
  reldistMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "reldist.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools reldist"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void reldist_help(void);

int reldist_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;
    bool summary = true;

    // input arguments
    bool haveBedA           = false;
    bool haveBedB           = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) reldist_help();

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
        else if(PARAMETER_CHECK("-detail", 7, parameterLength)) {
            summary = false;
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
        cerr << endl 
             << "*****" 
             << endl 
             << "*****ERROR: Need -a and -b files. " 
             << endl 
             << "*****" 
             << endl;
        showHelp = true;
    }


    if (!showHelp) {

        RelativeDistance *rd = new RelativeDistance(bedAFile, 
                                                    bedBFile,
                                                    summary);
        delete rd;
        return 0;
    }
    else {
        reldist_help();
        return 0;
    }
}

void reldist_help(void) {

    cerr << "\nTool:    bedtools reldist" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Calculate the relative distance distribution "
         << "b/w two feature files."
         << endl << endl;

    cerr << "Usage:   " 
         << PROGRAM_NAME 
         << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" 
         << endl << endl;

    cerr << "Options: " << endl;


    cerr << "\t-detail\t"       << "Report the relative" 
                                << "distance for each interval in A" 
                                << endl << endl;
    // end the program here
    exit(1);

}
