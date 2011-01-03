/*****************************************************************************
  annotateMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "annotateBed.h"
#include "version.h"

using namespace std;

// define the version
#define PROGRAM_NAME "annotateBed"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input file
    string mainFile;

    // parm flags
    bool forceStrand    = false;
    bool haveBed        = false;
    bool haveFiles      = false;
    bool haveTitles     = false;
    bool reportCounts   = false;
    bool reportBoth     = false;

    // list of annotation files / names
    vector<string> inputFiles;
    vector<string> inputTitles;

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

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBed  = true;
                mainFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-files", 6, parameterLength)) {
            if ((i+1) < argc) {
                haveFiles = true;
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    inputFiles.push_back(file);
                    i++;
                    if (i < argc)
                        file = argv[i];
                }
                i--;
            }
        }
        else if(PARAMETER_CHECK("-names", 6, parameterLength)) {
            if ((i+1) < argc) {
                haveTitles = true;
                i = i+1;
                string title = argv[i];
                while (title[0] != '-' && i < argc) {
                    inputTitles.push_back(title);
                    i++;
                    if (i < argc)
                        title = argv[i];
                }
                i--;
            }
        }
        else if(PARAMETER_CHECK("-counts", 7, parameterLength)) {
            reportCounts = true;
        }
        else if(PARAMETER_CHECK("-both", 5, parameterLength)) {
            reportBoth = true;
        }
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            forceStrand = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBed || !haveFiles) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i and -files files. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedAnnotate *ba = new BedAnnotate(mainFile, inputFiles, inputTitles, forceStrand, reportCounts, reportBoth);
        ba->AnnotateBed();
        delete ba;
        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Annotates the depth & breadth of coverage of features from multiple files" << endl;
    cerr << "\t on the intervals in -i." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> -files FILE1 FILE2 .. FILEn" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-names\t"        << "A list of names (one / file) to describe each file in -i." << endl;
    cerr                        << "\t\tThese names will be printed as a header line." << endl << endl;

    cerr << "\t-counts\t"       << "Report the count of features in each file that overlap -i." << endl;
    cerr                        << "\t\t- Default is to report the fraction of -i covered by each file." << endl << endl;

    cerr << "\t-both\t"         << "Report the counts followed by the % coverage." << endl;
    cerr                        << "\t\t- Default is to report the fraction of -i covered by each file." << endl << endl;

    cerr << "\t-s\t"            << "Force strandedness.  That is, only include hits in A that" << endl;
    cerr                        << "\t\toverlap B on the same strand." << endl;
    cerr                        << "\t\t- By default, hits are included without respect to strand." << endl << endl;

    exit(1);
}
