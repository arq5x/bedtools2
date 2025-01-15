#include "revcompBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools revcomp"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void revcomp_help(int);

int revcomp_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile = "stdin";
    string genomeFile;

    bool haveBed    = true;
    bool haveGenome = false;
    bool printHeader = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) revcomp_help(0);

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

    if (showHelp) revcomp_help(1);

    BedRevcomp(bedFile, genomeFile, printHeader);
    return 0;
}

void revcomp_help(int exit_status) {

    cerr << "\nTool:    bedtools revcomp (aka revcompBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Compute the intervals on the reverse complement." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> -g <genome>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-header\t"           << "Print the header from the input file prior to results." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  Starts will be set to 0 if options would force it below 0." << endl;
    cerr << "\t(2)  Ends will be set to the chromosome length if above the max chrom length." << endl;
    cerr << "\t(3)  The genome file should be tab delimited and structured as follows:" << endl;
    cerr << "\n\t<chromName><TAB><chromSize>" << endl << endl;
    cerr << "\tFor example, Human (hg19):" << endl;
    cerr << "\tchr1\t249250621" << endl;
    cerr << "\tchr2\t243199373" << endl;
    cerr << "\t..." << endl;
    cerr << "\tchr18_gl000207_random\t4262" << endl << endl;

    cerr << "Tip 1. Use samtools faidx to create a genome file from a FASTA: " << endl;
    cerr << "\tOne can the samtools faidx command to index a FASTA file." << endl;
    cerr << "\tThe resulting .fai index is suitable as a genome file, " << endl;
    cerr << "\tas bedtools will only look at the first two, relevant columns" << endl;
    cerr << "\tof the .fai file." << endl << endl;
    cerr << "\tFor example:" << endl;
    cerr << "\tsamtools faidx GRCh38.fa" << endl;
    cerr << "\tbedtools revcomp -i my.bed -g GRCh38.fa.fai" << endl << endl;

    cerr << "Tip 2. Use UCSC Table Browser to create a genome file: " << endl;
    cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract" << endl;
    cerr << "\tchromosome sizes. For example, H. sapiens:" << endl << endl;
    cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \\" << endl;
    cerr << "\t\"select chrom, size from hg19.chromInfo\"  > hg19.genome" << endl << endl;

    // end the program here
    exit(exit_status);
}
