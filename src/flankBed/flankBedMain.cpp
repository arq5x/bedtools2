/*****************************************************************************
  flankBedMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "flankBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools flank"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void flank_help(void);

int flank_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile = "stdin";
    string genomeFile;

    bool haveBed    = true;
    bool haveGenome = false;
    bool haveLeft   = false;
    bool haveRight  = false;
    bool haveBoth   = false;

    bool forceStrand = false;
    float leftSlop   = 0.0;
    float rightSlop  = 0.0;
    bool  fractional = false;
    bool printHeader        = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) flank_help();

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
        else if(PARAMETER_CHECK("-l", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveLeft = true;
                leftSlop = atof(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-r", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveRight = true;
                rightSlop = atof(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBoth = true;
                leftSlop = atof(argv[i + 1]);
                rightSlop = atof(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
            forceStrand = true;
        }
        else if(PARAMETER_CHECK("-pct", 4, parameterLength)) {
            fractional = true;
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
    if (!haveBed || !haveGenome) {
      cerr << endl << "*****" << endl << "*****ERROR: Need both a BED (-i) and a genome (-g) file. " << endl << "*****" << endl;
      showHelp = true;
    }
    if (!haveLeft && !haveRight && !haveBoth) {
      cerr << endl << "*****" << endl << "*****ERROR: Need -l and -r together or -b alone. " << endl << "*****" << endl;
      showHelp = true;
    }
    if ((!haveLeft && haveRight) || (haveLeft && !haveRight)) {
      cerr << endl << "*****" << endl << "*****ERROR: Need both -l and -r. " << endl << "*****" << endl;
      showHelp = true;
    }
    if (forceStrand && ((!(haveLeft) || !(haveRight)) && (!haveBoth))) {
      cerr << endl << "*****" << endl << "*****ERROR: Must supply -l and -r or just -b with -s. " << endl << "*****" << endl;
      showHelp = true;
    }

    if (!showHelp) {
        BedFlank *bc = new BedFlank(bedFile, genomeFile, forceStrand, 
                                    leftSlop, rightSlop, fractional, printHeader);
        delete bc;

        return 0;
    }
    else {
        flank_help();
    }
    return 0;
}

void flank_help(void) {
    
    cerr << "\nTool:    bedtools flank (aka flankBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Creates flanking interval(s) for each BED/GFF/VCF feature." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> -g <genome> [-b <int> or (-l and -r)]" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-b\t"                << "Create flanking interval(s) using -b base pairs in each direction." << endl;
    cerr                            << "\t\t- (Integer) or (Float, e.g. 0.1) if used with -pct." << endl << endl;

    cerr << "\t-l\t"                << "The number of base pairs that a flank should start from" << endl;
    cerr                            << "\t\torig. start coordinate." << endl;    
    cerr                            << "\t\t- (Integer) or (Float, e.g. 0.1) if used with -pct." << endl << endl;
        
    cerr << "\t-r\t"                << "The number of base pairs that a flank should end from" << endl;
    cerr                            << "\t\torig. end coordinate." << endl;
    cerr                            << "\t\t- (Integer) or (Float, e.g. 0.1) if used with -pct." << endl << endl;
        
    cerr << "\t-s\t"                << "Define -l and -r based on strand." << endl;
    cerr                            << "\t\tE.g. if used, -l 500 for a negative-stranded feature, " << endl;
    cerr                            << "\t\tit will start the flank 500 bp downstream.  Default = false." << endl << endl;

    cerr << "\t-pct\t"              << "Define -l and -r as a fraction of the feature's length." << endl;
    cerr                            << "\t\tE.g. if used on a 1000bp feature, -l 0.50, " << endl;
    cerr                            << "\t\twill add 500 bp \"upstream\".  Default = false." << endl << endl;
    
    cerr << "\t-header\t"           << "Print the header from the input file prior to results." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  Starts will be set to 0 if options would force it below 0." << endl;
    cerr << "\t(2)  Ends will be set to the chromosome length if requested flank would" << endl;
    cerr <<        "\tforce it above the max chrom length." << endl;
    cerr << "\t(3)  In contrast to slop, which _extends_ intervals, bedtools flank" << endl;
    cerr <<        "\tcreates new intervals from the regions just up- and down-stream" << endl;
    cerr <<        "\tof your existing intervals." << endl;

    cerr << "\t(4)  The genome file should tab delimited and structured as follows:" << endl;
    cerr << "\n\t<chromName><TAB><chromSize>" << endl << endl;
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
