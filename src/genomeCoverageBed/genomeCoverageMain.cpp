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
#define PROGRAM_NAME "bedtools genomecov"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void genomecoverage_help(void);

int genomecoverage_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile;
    string genomeFile;
    int max = INT_MAX;
    float scale = 1.0;

    bool haveBed = false;
    bool bamInput = false;
    bool haveGenome = false;
    bool startSites = false;
    bool bedGraph = false;
    bool bedGraphAll = false;
    bool eachBase = false;
    bool eachBaseZeroBased = false;
    bool obeySplits = false;
    bool haveScale = false;
    bool filterByStrand = false;
    bool only_5p_end = false;
    bool only_3p_end = false;
    bool add_gb_track_line = false;
    string gb_track_opts;
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

    if(showHelp) genomecoverage_help();

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
                haveBed = true;
                bamInput = true;
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
        else if(PARAMETER_CHECK("-d", 2, parameterLength)) {
            eachBase = true;
        }
        else if(PARAMETER_CHECK("-dz", 3, parameterLength)) {
            eachBase = true;
            eachBaseZeroBased = true;
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
        else if(PARAMETER_CHECK("-scale", 6, parameterLength)) {
            if ((i+1) < argc) {
                haveScale = true;
                scale = atof(argv[i + 1]);
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
        else if(PARAMETER_CHECK("-3", 2, parameterLength)) {
                only_3p_end = true;
        }
        else if(PARAMETER_CHECK("-5", 2, parameterLength)) {
                only_5p_end = true;
        }
        else if(PARAMETER_CHECK("-trackline", 10, parameterLength)) {
                add_gb_track_line = true;
        }
        else if(PARAMETER_CHECK("-trackopts", 10, parameterLength)) {
                if ((i+1) < argc) {
                    add_gb_track_line = true;
                    gb_track_opts = argv[i+1];
                    i++;
                } else {
                    cerr << "*****ERROR: -trackopts options requires a value (UCSC/GB track definition parameters)" << endl;
                    showHelp = true;
                }
        }
        else {
          cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBed && !haveGenome && !bamInput) {
      cerr << endl << "*****" << endl << "*****ERROR: Need both a BED (-i) and a genome (-g) file. " << endl << "*****" << endl;
      showHelp = true;
    }
    if (bedGraph && eachBase) {
      cerr << endl << "*****" << endl << "*****ERROR: Use -d/-dz or -bg, not both" << endl << "*****" << endl;
      showHelp = true;
    }
    if (bedGraphAll && eachBase) {
      cerr << endl << "*****" << endl << "*****ERROR: Use -d/-dz or -bga, not both" << endl << "*****" << endl;
      showHelp = true;
    }

    if (only_3p_end && only_5p_end) {
      cerr << endl << "*****" << endl << "*****ERROR: Use -3 or -5, not both " << endl << "*****" << endl;
      showHelp = true;
    }

    if ( (only_3p_end||only_5p_end) && obeySplits) {
      cerr << endl << "*****" << endl << "*****ERROR: Use -split can't be used with -3 or -5." << endl << "*****" << endl;
      showHelp = true;
    }

    if (add_gb_track_line && !(bedGraph||bedGraphAll)) {
      cerr << endl << "*****" << endl << "*****ERROR: Using -trackline requires bedGraph output (use -bg or -bga)." << endl << "*****" << endl;
      showHelp = true;
    }

    if (haveScale && !(bedGraph||bedGraphAll||eachBase)) {
      cerr << endl << "*****" << endl << "*****ERROR: Using -scale requires bedGraph output (use -bg or -bga) or per base depth (-d)." << endl << "*****" << endl;
      showHelp = true;
    }
    
    if (!showHelp) {
        BedGenomeCoverage *bc = new BedGenomeCoverage(bedFile, genomeFile, eachBase,
                                                      startSites, bedGraph, bedGraphAll,
                                                      max, scale, bamInput, obeySplits,
                                                      filterByStrand, requestedStrand,
                                                      only_5p_end, only_3p_end,
                                                      eachBaseZeroBased,
                                                      add_gb_track_line, gb_track_opts);
        delete bc;
    }
    else {
        genomecoverage_help();
    }
    return 0;
}

void genomecoverage_help(void) {

    cerr << "\nTool:    bedtools genomecov (aka genomeCoverageBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Compute the coverage of a feature file among a genome." << endl << endl;

    cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> -g <genome>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-ibam\t\t" << "The input file is in BAM format." << endl;
    cerr << "\t\t\tNote: BAM _must_ be sorted by position" << endl << endl;

    cerr << "\t-d\t\t" << "Report the depth at each genome position (with one-based coordinates)." << endl;
    cerr << "\t\t\tDefault behavior is to report a histogram." << endl << endl;

    cerr << "\t-dz\t\t" << "Report the depth at each genome position (with zero-based coordinates)." << endl;
    cerr << "\t\t\tReports only non-zero positions." << endl;
    cerr << "\t\t\tDefault behavior is to report a histogram." << endl << endl;

    cerr << "\t-bg\t\t" << "Report depth in BedGraph format. For details, see:" << endl;
    cerr << "\t\t\tgenome.ucsc.edu/goldenPath/help/bedgraph.html" << endl << endl;

    cerr << "\t-bga\t\t" << "Report depth in BedGraph format, as above (-bg)." << endl;
    cerr << "\t\t\tHowever with this option, regions with zero " << endl;
    cerr << "\t\t\tcoverage are also reported. This allows one to" << endl;
    cerr << "\t\t\tquickly extract all regions of a genome with 0 " << endl;
    cerr << "\t\t\tcoverage by applying: \"grep -w 0$\" to the output." << endl << endl;

    cerr << "\t-split\t\t" << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl;
    cerr << "\t\t\twhen computing coverage." << endl;
    cerr << "\t\t\tFor BAM files, this uses the CIGAR \"N\" and \"D\" operations " << endl;
    cerr << "\t\t\tto infer the blocks for computing coverage." << endl;
    cerr << "\t\t\tFor BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds" << endl;
    cerr << "\t\t\tfields (i.e., columns 10,11,12)." << endl << endl;

    cerr << "\t-strand\t\t" << "Calculate coverage of intervals from a specific strand." << endl;
    cerr << "\t\t\tWith BED files, requires at least 6 columns (strand is column 6). " << endl;
    cerr << "\t\t\t- (STRING): can be + or -" << endl << endl;

    cerr << "\t-5\t\t" << "Calculate coverage of 5\" positions (instead of entire interval)." << endl << endl;

    cerr << "\t-3\t\t" << "Calculate coverage of 3\" positions (instead of entire interval)." << endl << endl;

    cerr << "\t-max\t\t" << "Combine all positions with a depth >= max into" << endl;
    cerr << "\t\t\ta single bin in the histogram. Irrelevant" << endl;
    cerr << "\t\t\tfor -d and -bedGraph" << endl;
    cerr << "\t\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-scale\t\t" << "Scale the coverage by a constant factor." << endl;
    cerr << "\t\t\tEach coverage value is multiplied by this factor before being reported." << endl;
    cerr << "\t\t\tUseful for normalizing coverage by, e.g., reads per million (RPM)." << endl;
    cerr << "\t\t\t- Default is 1.0; i.e., unscaled." << endl;
    cerr << "\t\t\t- (FLOAT)" << endl << endl;

    cerr << "\t-trackline\t" << "Adds a UCSC/Genome-Browser track line definition in the first line of the output." << endl;
    cerr <<"\t\t\t- See here for more details about track line definition:" << endl;
    cerr <<"\t\t\t      http://genome.ucsc.edu/goldenPath/help/bedgraph.html" << endl;
    cerr <<"\t\t\t- NOTE: When adding a trackline definition, the output BedGraph can be easily" << endl;
    cerr <<"\t\t\t      uploaded to the Genome Browser as a custom track," << endl;
    cerr <<"\t\t\t      BUT CAN NOT be converted into a BigWig file (w/o removing the first line)." << endl << endl;

    cerr << "\t-trackopts\t"<<"Writes additional track line definition parameters in the first line." << endl;
    cerr <<"\t\t\t- Example:" << endl;
    cerr <<"\t\t\t   -trackopts 'name=\"My Track\" visibility=2 color=255,30,30'" << endl;
    cerr <<"\t\t\t   Note the use of single-quotes if you have spaces in your parameters." << endl;
    cerr <<"\t\t\t- (TEXT)" << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1) The genome file should tab delimited and structured as follows:" << endl;
    cerr << "\t <chromName><TAB><chromSize>" << endl << endl;
    cerr << "\tFor example, Human (hg19):" << endl;
    cerr << "\tchr1\t249250621" << endl;
    cerr << "\tchr2\t243199373" << endl;
    cerr << "\t..." << endl;
    cerr << "\tchr18_gl000207_random\t4262" << endl << endl;

    cerr << "\t(2) The input BED (-i) file must be grouped by chromosome." << endl;
    cerr << "\t A simple \"sort -k 1,1 <BED> > <BED>.sorted\" will suffice."<< endl << endl;

    cerr << "\t(3) The input BAM (-ibam) file must be sorted by position." << endl;
    cerr << "\t A \"samtools sort <BAM>\" should suffice."<< endl << endl;

    cerr << "Tips: " << endl;
    cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract" << endl;
    cerr << "\tchromosome sizes. For example, H. sapiens:" << endl << endl;
    cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \\" << endl;
    cerr << "\t\"select chrom, size from hg19.chromInfo\" > hg19.genome" << endl << endl;


    // end the program here
    exit(1);
}

