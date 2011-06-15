/*****************************************************************************
  nucBedMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "nucBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "nucBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string fastaDbFile;
    string bedFile;
    string pattern;

    // checks for existence of parameters
    bool haveFastaDb = false;
    bool haveBed     = false;
    bool printSeq    = false;
    bool hasPattern  = false;

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

        if(PARAMETER_CHECK("-fi", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveFastaDb = true;
                fastaDbFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-bed", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveBed = true;
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-seq", 4, parameterLength)) {
            printSeq = true;
        }
        else if(PARAMETER_CHECK("-pattern", 8, parameterLength)) {
            if ((i+1) < argc) {
                hasPattern = true;
                pattern = argv[i + 1];
                i++;
            }
        }
        else {
            cerr << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    if (!haveFastaDb || !haveBed) {
        showHelp = true;
    }

    if (!showHelp) {

        NucBed *nuc = new NucBed(fastaDbFile, bedFile, printSeq, hasPattern, pattern);
        delete nuc;

        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Profiles the nucleotide content of intervals in a fasta file." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -fi <fasta> -bed <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-fi\tInput FASTA file" << endl << endl;
    cerr << "\t-bed\tBED/GFF/VCF file of ranges to extract from -fi" << endl << endl;
    cerr << "\t-seq\tPrint the extracted sequence" << endl << endl;
    cerr << "\t-pattern\tReport the number of times a user-defined sequence is observed (case-insensitive)." << endl << endl;
    
    
    cerr << "Output format: " << endl;
    cerr << "\tThe following information will be reported after each original BED entry:" << endl;
    cerr << "\t    1) %AT content" << endl;
    cerr << "\t    2) %GC content" << endl;
    cerr << "\t    3) Number of As observed" << endl;
    cerr << "\t    4) Number of Cs observed" << endl;
    cerr << "\t    5) Number of Gs observed" << endl;
    cerr << "\t    6) Number of Ts observed" << endl;
    cerr << "\t    7) Number of Ns observed" << endl;
    cerr << "\t    8) Number of other bases observed" << endl;
    cerr << "\t    9) The length of the explored sequence/interval." << endl;
    cerr << "\t    10) The sequence extracted from the FASTA file. (optional, if -seq is used)" << endl;
    cerr << "\t    11) The number of times a user defined pattern was observed. (optional, if -pattern is used.)" << endl;

    // end the program here
    exit(1);

}
