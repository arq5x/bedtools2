/*****************************************************************************
  annotateMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "tagBam.h"
#include "version.h"

using namespace std;

// define the version
#define PROGRAM_NAME "tagBam"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input file
    string bamFile;
    float overlapFraction = 1E-9;
    string tag = "YB";

    // parm flags
    bool haveTag        = false;
    bool haveFraction   = false;
    bool forceStrand    = false;
    bool haveBam        = false;
    bool haveFiles      = false;
    bool haveLabels     = false;

    // list of annotation files / names
    vector<string> inputFiles;
    vector<string> inputLabels;

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
                haveBam  = true;
                bamFile = argv[i + 1];
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
        else if(PARAMETER_CHECK("-labels", 7, parameterLength)) {
            if ((i+1) < argc) {
                haveLabels = true;
                i = i+1;
                string label = argv[i];
                while (label[0] != '-' && i < argc) {
                    inputLabels.push_back(label);
                    i++;
                    if (i < argc)
                        label = argv[i];
                }
                i--;
            }
        }
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            forceStrand = true;
        }
        else if(PARAMETER_CHECK("-f", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveFraction = true;
                overlapFraction = atof(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-tag", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveTag = true;
                tag = argv[i + 1];
                i++;
            }
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBam || !haveFiles || !haveLabels) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i, -files, and -labels. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (haveTag && tag.size() > 2) {
        cerr << endl << "*****" << endl << "*****ERROR: Custom tags should be at most two characters per the SAM specification. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        TagBam *ba = new TagBam(bamFile, inputFiles, inputLabels, tag, forceStrand, overlapFraction);
        ba->Tag();
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

    cerr << "Summary: Annotates a BAM file based on overlaps with multiple BED/GFF/VCF files" << endl;
    cerr << "\t on the intervals in -i." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <BAM> -files FILE1 .. FILEn  -labels LAB1 .. LABn" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-s\t"            << "Force strandedness.  That is, only tag alignments that have the same" << endl;
    cerr                        << "\t\tstrand as a feature in the annotation file(s)." << endl << endl;

    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of the alignment." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-tag\t"          << "Dictate what the tag should be. Default is YB." << endl;
    cerr                        << "\t\t- STRING (two characters, e.g., YK)" << endl << endl;
    
    exit(1);
}
