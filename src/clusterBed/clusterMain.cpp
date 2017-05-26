/*****************************************************************************
  clusterMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "clusterBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools cluster"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void cluster_help(void);

int cluster_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile;
    int maxDistance = 0;

    // input arguments
    bool haveBed         = false;
    bool haveMaxDistance = false;
    bool forceStrand     = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) cluster_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedFile = argv[i + 1];
                i++;
                haveBed = true;
            }
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
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have input 
    if (!haveBed) {
        if (!isatty(STDIN_FILENO))
        {
            bedFile = "stdin";
        }
        else 
        {
            cerr << endl << "*****" << endl << "*****ERROR: Need -i BED file. " << endl << "*****" << endl;
            showHelp = true;
        }
    }

    if (!showHelp) {
        BedCluster *bc = new BedCluster(bedFile, maxDistance, forceStrand);
        delete bc;
    }
    else {
        cluster_help();
    }
    return 0;
}

void cluster_help(void) {
    
    cerr << "\nTool:    bedtools cluster" << endl;
    cerr << "Version: " << VERSION << "\n";        
    cerr << "Summary: Clusters overlapping/nearby BED/GFF/VCF intervals." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-s\t"                     << "Force strandedness.  That is, only merge features" << endl;
    cerr                                 << "\t\tthat are the same strand." << endl;
    cerr                                 << "\t\t- By default, merging is done without respect to strand." << endl << endl;

    cerr << "\t-d\t"                     << "Maximum distance between features allowed for features" << endl;
    cerr                                 << "\t\tto be merged." << endl;
    cerr                                 << "\t\t- Def. 0. That is, overlapping & book-ended features are merged." << endl;
    cerr                                 << "\t\t- (INTEGER)" << endl << endl;
    
    // end the program here
    exit(1);

}
