/*****************************************************************************
  randomBedMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "randomBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools random"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void random_help(void);

int random_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile = "stdin";
    string genomeFile;

    bool haveGenome       = false;
    bool haveSeed         = false;
    int seed              = -1;
    int length            = 100;
    int numToGenerate     = 1000000;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) random_help();

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
        else if(PARAMETER_CHECK("-seed", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveSeed = true;
                seed = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-l", 2, parameterLength)) {
            if ((i+1) < argc) {
                length = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-n", 2, parameterLength)) {
            if ((i+1) < argc) {
                numToGenerate = atoi(argv[i + 1]);
                i++;
            }
        }
        else {
          cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveGenome) {
      cerr << endl << "*****" << endl << "*****ERROR: Need a genome (-g) file. " << endl << "*****" << endl;
      showHelp = true;
    }

    if (!showHelp) {
        BedRandom *br = new BedRandom(genomeFile, numToGenerate, seed, haveSeed, length);
        delete br;
        return 0;
    }
    else {
        random_help();
    }
    return 0;
}

void random_help(void) {

    cerr << "\nTool:    bedtools random (aka randomBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Generate random intervals among a genome." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -g <genome>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-l\t"                << "The length of the intervals to generate." << endl;
    cerr                            << "\t\t- Default = 100." << endl;
    cerr                            << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-n\t"                << "The number of intervals to generate." << endl;
    cerr                            << "\t\t- Default = 1,000,000." << endl;
    cerr                            << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-seed\t"             << "Supply an integer seed for the shuffling." << endl;
    cerr                            << "\t\t- By default, the seed is chosen automatically." << endl;
    cerr                            << "\t\t- (INTEGER)" << endl << endl;

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


    // end the program here
    exit(1);
}
