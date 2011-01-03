/*****************************************************************************
  genomeCoverageMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "genomeCoverageBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "genomeCoverageBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile;
    string genomeFile;
    int max = INT_MAX;

    bool haveBed           = false;
    bool bamInput          = false;
    bool haveGenome        = false;
    bool startSites        = false;
    bool bedGraph          = false;
    bool bedGraphAll       = false;
    bool eachBase          = false;
    bool obeySplits        = false;
    bool filterByStrand    = false;
    string requestedStrand = "X";

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
                haveBed = true;
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-ibam", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveBed  = true;
                bamInput = true;
                bedFile  = argv[i + 1];
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
        else if(PARAMETER_CHECK("-d", 2, parameterLength)) {
            eachBase = true;
        }
        else if(PARAMETER_CHECK("-bg", 3, parameterLength)) {
            bedGraph = true;
        }
        else if(PARAMETER_CHECK("-bga", 4, parameterLength)) {
            bedGraphAll = true;
        }
        else if(PARAMETER_CHECK("-max", 4, parameterLength)) {
            if ((i+1) < argc) {
                max = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-split", 6, parameterLength)) {
            obeySplits = true;
        }
        else if(PARAMETER_CHECK("-strand", 7, parameterLength)) {
            if ((i+1) < argc) {
                filterByStrand = true;
                requestedStrand = argv[i+1][0];
                if (!(requestedStrand == "-" || requestedStrand == "+")) {
                    cerr << "*****ERROR: invalid -strand value (" << requestedStrand << "). Allowed options are + or -" << endl;
                    showHelp = true;
                }
                i++;
            }
            else {
                cerr << "*****ERROR: -strand options requires a value: + or -" << endl;
                showHelp = true;
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
    if (bedGraph && eachBase) {
      cerr << endl << "*****" << endl << "*****ERROR: Use -d or -bg, not both" << endl << "*****" << endl;
      showHelp = true;
    }
    if (bedGraphAll && eachBase) {
      cerr << endl << "*****" << endl << "*****ERROR: Use -d or -bga, not both" << endl << "*****" << endl;
      showHelp = true;
    }

    if (!showHelp) {

        BedGenomeCoverage *bc = new BedGenomeCoverage(bedFile, genomeFile, eachBase,
                                                      startSites, bedGraph, bedGraphAll,
                                                      max, bamInput, obeySplits,
                                                      filterByStrand, requestedStrand);
        delete bc;

        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Authors: Aaron Quinlan (aaronquinlan@gmail.com)" << endl;
    cerr << "         Assaf Gordon, CSHL" << endl << endl;

    cerr << "Summary: Compute the coverage of a feature file among a genome." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> -g <genome>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-ibam\t"         << "The input file is in BAM format." << endl;
    cerr                        << "\t\tNote: BAM _must_ be sorted by position" << endl << endl;

    cerr << "\t-d\t"            << "Report the depth at each genome position." << endl;
    cerr                        << "\t\tDefault behavior is to report a histogram." << endl << endl;

    cerr << "\t-bg\t"           << "Report depth in BedGraph format. For details, see:" << endl;
    cerr                        << "\t\tgenome.ucsc.edu/goldenPath/help/bedgraph.html" << endl << endl;

    cerr << "\t-bga\t"          << "Report depth in BedGraph format, as above (-bg)." << endl;
    cerr                        << "\t\tHowever with this option, regions with zero " << endl;
    cerr                        << "\t\tcoverage are also reported.  This allows one to" << endl;
    cerr                        << "\t\tquickly extract all regions of a genome with 0 " << endl;
    cerr                        << "\t\tcoverage by applying: \"grep -w 0$\" to the output." << endl << endl;

    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl;
    cerr                        << "\t\twhen computing coverage." << endl;
    cerr                        << "\t\tFor BAM files, this uses the CIGAR \"N\" and \"D\" operations " << endl;
    cerr                        << "\t\tto infer the blocks for computing coverage." << endl;
    cerr                        << "\t\tFor BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds" << endl;
    cerr                        << "\t\tfields (i.e., columns 10,11,12)." << endl << endl;

    cerr << "\t-strand\t"       << "Calculate coverage of intervals from a specific strand." << endl;
    cerr                        << "\t\tWith BED files, requires at least 6 columns (strand is column 6). " << endl;
    cerr                        << "\t\t- (STRING): can be + or -" << endl << endl;


    cerr << "\t-max\t"          << "Combine all positions with a depth >= max into" << endl;
    cerr                        << "\t\ta single bin in the histogram. Irrelevant" << endl;
    cerr                        << "\t\tfor -d and -bedGraph" << endl;
    cerr                        << "\t\t- (INTEGER)" << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  The genome file should tab delimited and structured as follows:" << endl;
    cerr << "\t     <chromName><TAB><chromSize>" << endl << endl;
    cerr << "\tFor example, Human (hg19):" << endl;
    cerr << "\tchr1\t249250621" << endl;
    cerr << "\tchr2\t243199373" << endl;
    cerr << "\t..." << endl;
    cerr << "\tchr18_gl000207_random\t4262" << endl << endl;

    cerr << "\t(2)  The input BED (-i) file must be grouped by chromosome." << endl;
    cerr << "\t     A simple \"sort -k 1,1 <BED> > <BED>.sorted\" will suffice."<< endl << endl;

    cerr << "\t(3)  The input BAM (-ibam) file must be sorted by position." << endl;
    cerr << "\t     A \"samtools sort <BAM>\" should suffice."<< endl << endl;

    cerr << "Tips: " << endl;
    cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract" << endl;
    cerr << "\tchromosome sizes. For example, H. sapiens:" << endl << endl;
    cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e /" << endl;
    cerr << "\t\"select chrom, size from hg19.chromInfo\" > hg19.genome" << endl << endl;


    // end the program here
    exit(1);
}
