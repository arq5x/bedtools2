/*****************************************************************************
  shuffleBedMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "shuffleBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools shuffle"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void shuffle_help(void);

int shuffle_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile = "stdin";
    string excludeFile;
    string includeFile;
    string genomeFile;

    bool haveBed          = true;
    bool haveGenome       = false;
    bool haveExclude      = false;
    bool haveInclude      = false;
    bool haveSeed         = false;
    float overlapFraction = 0.0;
    int seed              = -1;
    bool sameChrom        = false;
    bool chooseChrom      = false;


    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) shuffle_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-g", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveGenome = true;
                genomeFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-excl", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveExclude = true;
                excludeFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-incl", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveInclude = true;
                chooseChrom = true;
                includeFile = argv[i + 1];
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
        else if(PARAMETER_CHECK("-chrom", 6, parameterLength)) {
            chooseChrom = true;
            sameChrom = true;
        }
        else if(PARAMETER_CHECK("-chromFirst", 11, parameterLength)) {
            chooseChrom = true;
        }
        else if(PARAMETER_CHECK("-f", 2, parameterLength)) {
            if ((i+1) < argc) {
                overlapFraction = atof(argv[i + 1]);
                i++;
            }
        }
        else {
          cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBed || !haveGenome) {
      cerr << endl << "*****" << endl << "*****ERROR: Need both a BED (-i) and a genome (-g) file. " << endl << "*****" << endl;
      showHelp = true;
    }
    
    if (haveInclude && haveExclude) {
      cerr << endl << "*****" << endl << "*****ERROR: Cannot use -incl and -excl together." << endl << "*****" << endl;
      showHelp = true;
    }

    if (!showHelp) {
        BedShuffle *bc = new BedShuffle(bedFile, genomeFile, excludeFile, includeFile, 
                                        haveSeed, haveExclude, haveInclude, sameChrom, 
                                        overlapFraction, seed, chooseChrom);
        delete bc;
        return 0;
    }
    else {
        shuffle_help();
    }
    return 0;
}

void shuffle_help(void) {

    cerr << "\nTool:    bedtools shuffle (aka shuffleBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Randomly permute the locations of a feature file among a genome." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> -g <genome>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-excl\t"             << "A BED/GFF/VCF file of coordinates in which features in -i" << endl;
    cerr                            << "\t\tshould not be placed (e.g. gaps.bed)." << endl << endl;

    cerr << "\t-incl\t"             << "Instead of randomly placing features in a genome, the -incl" << endl;
    cerr                            << "\t\toptions defines a BED/GFF/VCF file of coordinates in which " << endl;
    cerr                            << "\t\tfeatures in -i should be randomly placed (e.g. genes.bed). " << endl;
    cerr                            << "\t\t- NOTE: Forces use of -chromFirst (see below)." << endl << endl;
    
    cerr << "\t-chrom\t"            << "Keep features in -i on the same chromosome."<< endl;
    cerr                            << "\t\t- By default, the chrom and position are randomly chosen." << endl;
    cerr                            << "\t\t- NOTE: Forces use of -chromFirst (see below)." << endl << endl;

    cerr << "\t-seed\t"             << "Supply an integer seed for the shuffling." << endl;
    cerr                            << "\t\t- By default, the seed is chosen automatically." << endl;
    cerr                            << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-f\t"                << "Maximum overlap (as a fraction of the -i feature) with an -excl" << endl;
    cerr                            << "\t\tfeature that is tolerated before searching for a new, " << endl;
    cerr                            << "\t\trandomized locus. For example, -f 0.10 allows up to 10%" << endl;
    cerr                            << "\t\tof a randomized feature to overlap with a given feature" << endl;
    cerr                            << "\t\tin the -excl file. **Cannot be used with -incl file.**" << endl;
    cerr                            << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                            << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-chromFirst\t"       << "\n\t\tInstead of choosing a position randomly among the entire" << endl;
    cerr                            << "\t\tgenome (the default), first choose a chrom randomly, and then" << endl;
    cerr                            << "\t\tchoose a random start coordinate on that chrom.  This leads" << endl; 
    cerr                            << "\t\tto features being ~uniformly distributed among the chroms," << endl; 
    cerr                            << "\t\tas opposed to features being distribute as a function of chrom size." << endl << endl; 


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
