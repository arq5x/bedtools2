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
#define PROGRAM_NAME "bedtools multicov"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void multibamcov_help(void);

int multibamcov_main(int argc, char* argv[]) {

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
    bool obeySplits        = false;
    bool sameStrand        = false;
    bool diffStrand        = false;
    float overlapFraction = 1E-9;
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

    if(showHelp) multibamcov_help();

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
        else if(PARAMETER_CHECK("-split", 6, parameterLength)) {
            obeySplits = true;
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
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            sameStrand = true;
        }
        else if (PARAMETER_CHECK("-S", 2, parameterLength)) {
            diffStrand = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << 
                argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    if (!showHelp) {
        MultiCovBam *mc = new MultiCovBam(bamFiles, bedFile, 
                                          minQual, properOnly, 
                                          keepDuplicates, keepFailedQC,
                                          obeySplits, sameStrand,
                                          diffStrand, overlapFraction,
                                          reciprocalFraction);
        mc->CollectCoverage();
        delete mc;
    }
    else {
        multibamcov_help();
    }
    return 0;
}

void multibamcov_help(void) {

    cerr << "\nTool:    bedtools multicov (aka multiBamCov)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Counts sequence coverage for multiple bams at specific loci." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -bams aln.1.bam aln.2.bam ... aln.n.bam -bed <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-bams\t"        << "The bam files." << endl << endl;

    cerr << "\t-bed\t"         << "The bed file." << endl << endl;
    
    cerr << "\t-split\t"       << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl << endl;
                               
    cerr << "\t-s\t"           << "Require same strandedness.  That is, only report hits in B" << endl;
    cerr                       << "\t\tthat overlap A on the _same_ strand." << endl;
    cerr                       << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;
                               
    cerr << "\t-S\t"           << "Require different strandedness.  That is, only report hits in B" << endl;
    cerr                       << "\t\tthat overlap A on the _opposite_ strand." << endl;
    cerr                       << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of each -bed record." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-r\t"            << "Require that the fraction overlap be reciprocal for each -bed and B." << endl;
    cerr                        << "\t\t- In other words, if -f is 0.90 and -r is used, this requires" << endl;
    cerr                        << "\t\t  that B overlap 90% of each -bed and the -bed record _also_ overlaps 90% of B." << endl << endl;

    cerr << "\t-q\t"           << "Minimum mapping quality allowed. Default is 0." << endl << endl;

    cerr << "\t-D\t"           << "Include duplicate reads.  Default counts non-duplicates only" << endl << endl;

    cerr << "\t-F\t"           << "Include failed-QC reads.  Default counts pass-QC reads only" << endl << endl;

    cerr << "\t-p\t"           << "Only count proper pairs.  Default counts all alignments with" << endl;
    cerr << "\t\t"             << "MAPQ > -q argument, regardless of the BAM FLAG field." << endl << endl;

    // end the program here
    exit(1);

}
