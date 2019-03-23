/*****************************************************************************
  windowMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "windowBed.h"
#include "version.h"

using namespace std;

// define the version
#define PROGRAM_NAME "bedtools window"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void window_help(void);


int window_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;

    // input arguments
    int leftSlop  = 1000;
    int rightSlop = 1000;

    bool haveBedA            = false;
    bool haveBedB            = false;
    bool noHit               = false;
    bool anyHit              = false;
    bool writeCount          = false;
    bool haveSlop            = false;
    bool haveLeft            = false;
    bool haveRight           = false;
    bool strandWindows       = false;
    bool matchOnSameStrand   = false;
    bool matchOnDiffStrand   = false;
    bool inputIsBam          = false;
    bool outputIsBam         = true;
    bool uncompressedBam     = false;
    bool printHeader         = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) window_help();

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
        else if(PARAMETER_CHECK("-abam", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveBedA = true;
                inputIsBam = true;
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
        else if(PARAMETER_CHECK("-bed", 4, parameterLength)) {
            outputIsBam = false;
        }
        else if(PARAMETER_CHECK("-u", 2, parameterLength)) {
            anyHit = true;
        }
        else if(PARAMETER_CHECK("-c", 2, parameterLength)) {
            writeCount = true;
        }
        else if (PARAMETER_CHECK("-v", 2, parameterLength)) {
            noHit = true;
        }
        else if (PARAMETER_CHECK("-sw", 3, parameterLength)) {
            strandWindows = true;
        }
        else if (PARAMETER_CHECK("-sm", 3, parameterLength)) {
            matchOnSameStrand = true;
        }
        else if (PARAMETER_CHECK("-Sm", 3, parameterLength)) {
            matchOnDiffStrand = true;
        }
        else if (PARAMETER_CHECK("-w", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveSlop = true;
                leftSlop = atoi(argv[i + 1]);
                rightSlop = leftSlop;
                i++;
            }
        }
        else if (PARAMETER_CHECK("-l", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveLeft = true;
                leftSlop = atoi(argv[i + 1]);
                i++;
            }
        }
        else if (PARAMETER_CHECK("-r", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveRight = true;
                rightSlop = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-ubam", 5, parameterLength)) {
            uncompressedBam = true;
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

    if (anyHit && noHit) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -v, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (anyHit && writeCount) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -c, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (haveLeft && (leftSlop < 0)) {
        cerr << endl << "*****" << endl << "*****ERROR: Upstream window (-l) must be positive." << endl << "*****" << endl;
        showHelp = true;
    }

    if (haveRight && (rightSlop < 0)) {
        cerr << endl << "*****" << endl << "*****ERROR: Downstream window (-r) must be positive." << endl << "*****" << endl;
        showHelp = true;
    }

    if (haveSlop && (haveLeft || haveRight)) {
        cerr << endl << "*****" << endl << "*****ERROR: Cannot choose -w with -l or -r.  Either specify -l and -r or specify solely -w" << endl << "*****" << endl;
        showHelp = true;
    }

    if ((haveLeft && !haveRight) || (haveRight && !haveLeft)) {
        cerr << endl << "*****" << endl << "*****ERROR: Please specify both -l and -r." << endl << "*****" << endl;
        showHelp = true;
    }

    if (matchOnSameStrand && matchOnDiffStrand) {
        cerr << endl << "*****" << endl << "*****ERROR: Use either -sm or -Sm, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedWindow *bi = new BedWindow(bedAFile, bedBFile, leftSlop, rightSlop, anyHit,
                                      noHit, writeCount, strandWindows, matchOnSameStrand, matchOnDiffStrand,
                                      inputIsBam, outputIsBam, uncompressedBam, printHeader);
        delete bi;
        return 0;
    }
    else {
        window_help();
    }
    return 0;
}


void window_help(void) {

    cerr << "\nTool:    bedtools window (aka windowBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Examines a \"window\" around each feature in A and" << endl;
    cerr << "\t reports all features in B that overlap the window. For each" << endl;
    cerr << "\t overlap the entire entry in A and B are reported." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-abam\t"         << "The A input file is in BAM format.  Output will be BAM as well. Replaces -a." << endl << endl;

    cerr << "\t-ubam\t"         << "Write uncompressed BAM output. Default writes compressed BAM." << endl << endl;

    cerr << "\t-bed\t"          << "When using BAM input (-abam), write output as BED. The default" << endl;
    cerr                        << "\t\tis to write output in BAM when using -abam." << endl << endl;

    cerr << "\t-w\t"            << "Base pairs added upstream and downstream of each entry" << endl;
    cerr                        << "\t\tin A when searching for overlaps in B." << endl;
    cerr                        << "\t\t- Creates symmetrical \"windows\" around A." << endl;
    cerr                        << "\t\t- Default is 1000 bp." << endl;
    cerr                        << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-l\t"            << "Base pairs added upstream (left of) of each entry" << endl;
    cerr                        << "\t\tin A when searching for overlaps in B." << endl;
    cerr                        << "\t\t- Allows one to define asymmetrical \"windows\"." << endl;
    cerr                        << "\t\t- Default is 1000 bp." << endl;
    cerr                        << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-r\t"            << "Base pairs added downstream (right of) of each entry" << endl;
    cerr                        << "\t\tin A when searching for overlaps in B." << endl;
    cerr                        << "\t\t- Allows one to define asymmetrical \"windows\"." << endl;
    cerr                        << "\t\t- Default is 1000 bp." << endl;
    cerr                        << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-sw\t"           << "Define -l and -r based on strand.  For example if used, -l 500" << endl;
    cerr                        << "\t\tfor a negative-stranded feature will add 500 bp downstream." << endl;
    cerr                        << "\t\t- Default = disabled." << endl << endl;

    cerr << "\t-sm\t"           << "Only report hits in B that overlap A on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-Sm\t"           << "Only report hits in B that overlap A on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-u\t"            << "Write the original A entry _once_ if _any_ overlaps found in B." << endl;
    cerr                        << "\t\t- In other words, just report the fact >=1 hit was found." << endl << endl;

    cerr << "\t-c\t"            << "For each entry in A, report the number of overlaps with B." << endl;
    cerr                        << "\t\t- Reports 0 for A entries that have no overlap with B." << endl;
    cerr                        << "\t\t- Overlaps restricted by -w, -l, and -r." << endl << endl;

    cerr << "\t-v\t"            << "Only report those entries in A that have _no overlaps_ with B." << endl;
    cerr                        << "\t\t- Similar to \"grep -v.\"" << endl << endl;
    
    cerr << "\t-header\t"       << "Print the header from the A file prior to results." << endl << endl;

    // end the program here
    exit(1);
}
