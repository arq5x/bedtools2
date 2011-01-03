/*****************************************************************************
  maskFastaFromBedMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "maskFastaFromBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "maskFastaFromBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string fastaInFile;
    string bedFile;

    // output files
    string fastaOutFile;

    // checks for existence of parameters
    bool haveFastaIn = false;
    bool haveBed = false;
    bool haveFastaOut = false;
    bool softMask = false;

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
                haveFastaIn = true;
                fastaInFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-fo", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveFastaOut = true;
                fastaOutFile = argv[i + 1];
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
        else if(PARAMETER_CHECK("-soft", 5, parameterLength)) {
            softMask = true;
        }
        else {
            cerr << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    if (!haveFastaIn || !haveFastaOut || !haveBed) {
        showHelp = true;
    }

    if (!showHelp) {

        MaskFastaFromBed *maskFasta = new MaskFastaFromBed(fastaInFile, bedFile, fastaOutFile, softMask);
        delete maskFasta;
        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {



    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Mask a fasta file based on feature coordinates." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -fi <fasta> -out <fasta> -bed <bed/gff/vcf>" << endl << endl;

    cerr << "Options:" << endl;
    cerr << "\t-fi\tInput FASTA file" << endl;
    cerr << "\t-bed\tBED/GFF/VCF file of ranges to mask in -fi" << endl;
    cerr << "\t-fo\tOutput FASTA file" << endl;
    cerr << "\t-soft\tEnforce \"soft\" masking.  That is, instead of masking with Ns," << endl;
    cerr << "\t\tmask with lower-case bases." << endl;

    // end the program here
    exit(1);

}
