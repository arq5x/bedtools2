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
#define PROGRAM_NAME "bedtools tag"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void tagbam_help(void);

int tagbam_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input file
    string bamFile;
    float overlapFraction = 1E-9;
    string tag = "YB";

    // parm flags
    bool haveTag          = false;
    bool haveFraction     = false;
    bool useNames         = false;
    bool useScores        = false;
    bool useIntervals     = false;
    bool sameStrand       = false;
    bool diffStrand       = false;
    bool haveBam          = false;
    bool haveFiles        = false;
    bool haveLabels       = false;


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

    if(showHelp) tagbam_help();

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
        else if (PARAMETER_CHECK("-names", 6, parameterLength)) {
            useNames = true;
        }
        else if (PARAMETER_CHECK("-scores", 7, parameterLength)) {
            useScores = true;
        }
        else if (PARAMETER_CHECK("-intervals", 10, parameterLength)) {
            useIntervals = true;
        }
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            sameStrand = true;
        }
        else if (PARAMETER_CHECK("-S", 2, parameterLength)) {
            diffStrand = true;
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
    if (!haveBam || !haveFiles) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i, -files" << endl << "*****" << endl;
        showHelp = true;
    }
    if (!useNames && !haveLabels && !useScores) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -labels or -names or -scores" << endl << "*****" << endl;
        showHelp = true;
    }
    if (sameStrand && diffStrand) {
        cerr << endl << "*****" << endl << "*****ERROR: Use -s or -S, not both. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (haveLabels && useNames) {
        cerr << endl << "*****" << endl << "*****ERROR: Use -labels or -names, not both. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (useScores && useNames) {
        cerr << endl << "*****" << endl << "*****ERROR: Use -scores or -names, not both. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (useScores && useIntervals) {
        cerr << endl << "*****" << endl << "*****ERROR: Use -scores or -intervals, not both. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (useNames && useIntervals) {
        cerr << endl << "*****" << endl << "*****ERROR: Use -names or -intervals, not both. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (!haveLabels && useIntervals) {
        cerr << endl << "*****" << endl << "*****ERROR: Supply -labels when using -intervals. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (haveTag && tag.size() > 2) {
        cerr << endl << "*****" << endl << "*****ERROR: Custom tags should be at most two characters per the SAM specification. " << endl << "*****" << endl;
        showHelp = true;
    }


    if (!showHelp) {
        TagBam *ba = new TagBam(bamFile, inputFiles, inputLabels, 
                                tag, useNames, useScores,  
                                useIntervals, sameStrand, diffStrand, 
                                overlapFraction);
        ba->Tag();
        delete ba;
        return 0;
    }
    else {
        tagbam_help();
    }
    return 0;
}

void tagbam_help(void) {

    cerr << "\nTool:    bedtools tag (aka tagBam)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Annotates a BAM file based on overlaps with multiple BED/GFF/VCF files" << endl;
    cerr << "\t on the intervals in -i." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <BAM> -files FILE1 .. FILEn  -labels LAB1 .. LABn" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-s\t"            << "Require overlaps on the same strand.  That is, only tag alignments that have the same" << endl;
    cerr                        << "\t\tstrand as a feature in the annotation file(s)." << endl << endl;

    cerr << "\t-S\t"            << "Require overlaps on the opposite strand.  That is, only tag alignments that have the opposite" << endl;
    cerr                        << "\t\tstrand as a feature in the annotation file(s)." << endl << endl;

    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of the alignment." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-tag\t"          << "Dictate what the tag should be. Default is YB." << endl;
    cerr                        << "\t\t- STRING (two characters, e.g., YK)" << endl << endl;
    
    cerr << "\t-names\t"        << "Use the name field from the annotation files to populate tags." << endl;
    cerr                        << "\t\tBy default, the -labels values are used." << endl << endl;

    cerr << "\t-scores\t"       << "Use the score field from the annotation files to populate tags." << endl;
    cerr                        << "\t\tBy default, the -labels values are used." << endl << endl;    

    cerr << "\t-intervals\t"    << "Use the full interval (including name, score, and strand) to populate tags." << endl;
    cerr                        << "\t\t\tRequires the -labels option to identify from which file the interval came." << endl << endl;    
    
    exit(1);
}
