/*****************************************************************************
  windowMakerMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "windowMaker.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools makewindows"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void windowmaker_help(void);

int windowmaker_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string genomeFile;

    // parms
    uint32_t size = 0;
    uint32_t step = 0;
    
    bool haveGenome = false;
    bool haveSize   = false;
    
    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) windowmaker_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-g", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveGenome = true;
                genomeFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-w", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveSize = true;
                size = atoi(argv[i + 1]);
                step = size;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
            if ((i+1) < argc) {
                step = atoi(argv[i + 1]);
                i++;
            }
        }
        else {
          cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveGenome || !haveSize) {
      cerr << endl << "*****" << endl << "*****ERROR: Need -g (genome file) and -w (window size). " << endl << "*****" << endl;
      showHelp = true;
    }
    if (!showHelp) {
        WindowMaker *wm = new WindowMaker(genomeFile, size, step);
        delete wm;
    }
    else {
        windowmaker_help();
    }
    return 0;
}

void windowmaker_help(void) {

    cerr << "\nTool:    bedtools makewindows" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Makes adjacent and/or sliding windows across a genome." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -g <genome> -w <window_size>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-s\t"            << "Step size: i.e., how many base pairs to step before" << endl;
    cerr                        << "\t\tcreating a new window. Used to create \"sliding\" windows." << endl;
    cerr                        << "\t\t- Defaults to -w (non-sliding windows)." << endl << endl;
    
    cerr << "Notes: " << endl;
    cerr << "\t(1)  The genome file should tab delimited and structured as follows:" << endl;
    cerr << "\t     <chromName><TAB><chromSize>" << endl << endl;
    cerr << "\tFor example, Human (hg19):" << endl;
    cerr << "\tchr1\t249250621" << endl;
    cerr << "\tchr2\t243199373" << endl;
    cerr << "\t..." << endl;
    cerr << "\tchr18_gl000207_random\t4262" << endl << endl;

    cerr << "Tips: " << endl;
    cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract" << endl;
    cerr << "\tchromosome sizes. For example, H. sapiens:" << endl << endl;
    cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \\" << endl;
    cerr << "\t\"select chrom, size from hg19.chromInfo\"  > hg19.genome" << endl << endl;

    exit(1);

}
