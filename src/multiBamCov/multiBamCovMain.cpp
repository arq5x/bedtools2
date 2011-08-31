/*****************************************************************************
  multiBamCovMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "multiBamCov.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "multiBamCov"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile;
    vector<string> bamFiles;
    int minQual = 0;

    // input arguments
    bool haveBed           = false;
    bool haveBams          = false;
    bool properOnly        = false;
    bool keepDuplicates    = false;
    bool keepFailedQC      = false;
    
    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

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

        if(PARAMETER_CHECK("-bed", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveBed = true;
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-bams", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveBams = true;
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    bamFiles.push_back(file);
                    i++;
                    if (i < argc)
                        file = argv[i];
                }
                i--;
            }
        }
        else if(PARAMETER_CHECK("-q", 2, parameterLength)) {
            if ((i+1) < argc) {
                minQual = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-p", 2, parameterLength)) {
            properOnly = true;
        }
        else if(PARAMETER_CHECK("-D", 2, parameterLength)) {
            keepDuplicates = true;
        }
        
        else if(PARAMETER_CHECK("-F", 2, parameterLength)) {
            keepFailedQC = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    if (!showHelp) {
        MultiCovBam *mc = new MultiCovBam(bamFiles, bedFile, minQual, properOnly, keepDuplicates, keepFailedQC);
        mc->CollectCoverage();
        delete mc;
        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Counts sequence coverage for multiple bams at specific loci." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -bams aln.1.bam aln.2.bam ... aln.n.bam -bed <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-bams\t"        << "The bam files." << endl << endl;

    cerr << "\t-bed\t"         << "The bed file." << endl << endl;

    cerr << "\t-q\t"           << "Minimum mapping quality allowed. Default is 0." << endl << endl;

    cerr << "\t-D\t"           << "Include duplicate-marked reads.  Default is to count non-duplicates only" << endl << endl;

    cerr << "\t-F\t"           << "Include failed-QC reads.  Default is to count pass-QC reads only" << endl << endl;

    cerr << "\t-p\t"           << "Only count proper pairs.  Default is to count all alignments with MAPQ" << endl;
    cerr << "\t\t"             << "greater than the -q argument, regardless of the BAM FLAG field." << endl << endl;

    // end the program here
    exit(1);

}
