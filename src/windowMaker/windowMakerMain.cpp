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
    string inputFile;
    WindowMaker::INPUT_FILE_TYPE inputFileType = WindowMaker::GENOME_FILE;
    WindowMaker::ID_METHOD idMethod = WindowMaker::ID_NONE;

    // parms
    uint32_t size = 0;
    uint32_t step = 0;
    uint32_t count = 0;

    bool haveGenome = false;
    bool haveBed = false;
    bool haveStep = false;
    bool haveSize = false;
    bool haveCount = false;
    bool reverse = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 6, parameterLength))) {
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
                inputFile = argv[i + 1];
                inputFileType = WindowMaker::GENOME_FILE;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBed = true;
                inputFile = argv[i + 1];
                inputFileType = WindowMaker::BED_FILE;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-w", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveSize = true;
                size = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
            if ((i+1) < argc) {
                step = atoi(argv[i + 1]);
                haveStep = true;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-n", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveCount = true;
                count = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                if (strcmp(argv[i+1],"winnum")==0)
                    idMethod = WindowMaker::ID_WINDOW_NUMBER;
                else if (strcmp(argv[i+1],"srcwinnum")==0)
                    idMethod = WindowMaker::ID_SOURCE_ID_WINDOW_NUMBER;
                else if (strcmp(argv[i+1],"src")==0)
                    idMethod = WindowMaker::ID_SOURCE_ID;
                else {
                    cerr << endl << "*****ERROR: Invalid ID method (" << argv[i+1] << "). Possible values are: winnum, srcwinnum" << endl << endl ;
                    showHelp = true;
                }
                i++;
            }
        }
        else if(PARAMETER_CHECK("-reverse", 8, parameterLength)) {
            reverse = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveStep) {
        step = size;
    }
    if (!haveGenome && !haveBed) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -g (genome file) or -b (BED file) for interval source. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (haveGenome && haveBed) {
        cerr << endl << "*****" << endl << "*****ERROR: Can't combine -g (genome file) and -b (BED file). Please use one or the other." << endl << "*****" << endl;
        showHelp = true;
    }
    if (!haveSize && !haveCount) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -w (window size) or -n (number of windows). " << endl << "*****" << endl;
        showHelp = true;
    }
    if (haveSize && haveCount) {
        cerr << endl << "*****" << endl << "*****ERROR: Can't combine -w (window size) and -n (number of windows). Please use one or the other. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (haveStep && haveCount) {
        cerr << endl << "*****" << endl << "*****ERROR: Can't combine -s (step size) and -n (number of windows). Please use one or the other. " << endl << "*****" << endl;
        showHelp = true;
    }
    if (step <= 0 && !haveBed && !haveCount && haveGenome) {
        cerr << endl << "*****" << endl << "*****ERROR: The step (-s) option must be greater than zero. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        WindowMaker *wm = NULL;
        if (haveCount)
            wm = new WindowMaker(inputFile, idMethod,
                                 inputFileType, count, reverse);
        if (haveSize)
            wm = new WindowMaker(inputFile, idMethod,
                                 inputFileType, size, step, reverse);
        delete wm;
    }
    else {
        windowmaker_help();
    }
    return 0;
}

void windowmaker_help(void) {

    cerr << "\nTool: bedtools makewindows" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Makes adjacent or sliding windows across a genome or BED file." << endl << endl;

    cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] [-g <genome> OR -b <bed>]" << endl;
    cerr << " [ -w <window_size> OR -n <number of windows> ]" << endl << endl;

    cerr << "Input Options: " << endl;

    cerr << "\t-g <genome>" << endl;
    cerr << "\t\tGenome file size (see notes below)." << endl;
    cerr << "\t\tWindows will be created for each chromosome in the file." << endl << endl;

    cerr << "\t-b <bed>" << endl;
    cerr << "\t\tBED file (with chrom,start,end fields)." << endl;
    cerr << "\t\tWindows will be created for each interval in the file." << endl << endl;

    cerr << "Windows Output Options: " << endl;

    cerr << "\t-w <window_size>" << endl;
    cerr << "\t\tDivide each input interval (either a chromosome or a BED interval)" << endl;
    cerr << "\t\tto fixed-sized windows (i.e. same number of nucleotide in each window)." << endl;
    cerr << "\t\tCan be combined with -s <step_size>" << endl << endl;

    cerr << "\t-s <step_size>" << endl;
    cerr << "\t\tStep size: i.e., how many base pairs to step before" << endl;
    cerr << "\t\tcreating a new window. Used to create \"sliding\" windows." << endl;
    cerr << "\t\t- Defaults to window size (non-sliding windows)." << endl << endl;

    cerr << "\t-n <number_of_windows>" << endl;
    cerr << "\t\tDivide each input interval (either a chromosome or a BED interval)" << endl;
    cerr << "\t\tto fixed number of windows (i.e. same number of windows, with" << endl;
    cerr << "\t\tvarying window sizes)." << endl << endl;

    cerr << "\t-reverse" << endl;
    cerr << "\t\t Reverse numbering of windows in the output, i.e. report " << endl;
    cerr << "\t\t windows in decreasing order" << endl << endl;

    cerr << "ID Naming Options: " << endl;
    cerr << "\t-i src|winnum|srcwinnum" << endl;
    cerr << "\t\tThe default output is 3 columns: chrom, start, end ." << endl;
    cerr << "\t\tWith this option, a name column will be added." << endl;
    cerr << "\t\t \"-i src\" - use the source interval's name." << endl;
    cerr << "\t\t \"-i winnum\" - use the window number as the ID (e.g. 1,2,3,4...)." << endl;
    cerr << "\t\t \"-i srcwinnum\" - use the source interval's name with the window number." << endl;
    cerr << "\t\tSee below for usage examples." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1) The genome file should tab delimited and structured as follows:" << endl;
    cerr << "\t <chromName><TAB><chromSize>" << endl << endl;
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
    cerr << "\tbedtools makewindows -w 100 -g GRCh38.fa.fai" << endl << endl;

    cerr << "Tip 2. Use UCSC Table Browser to create a genome file: " << endl;
    cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract" << endl;
    cerr << "\tchromosome sizes. For example, H. sapiens:" << endl << endl;
    cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \\" << endl;
    cerr << "\t\"select chrom, size from hg19.chromInfo\"  > hg19.genome" << endl << endl;

 
 
    cerr << "Examples: " << endl;
    cerr << " # Divide the human genome into windows of 1MB:" << endl;
    cerr << " $ " << PROGRAM_NAME << " -g hg19.txt -w 1000000" << endl;
    cerr << " chr1 0 1000000" << endl;
    cerr << " chr1 1000000 2000000" << endl;
    cerr << " chr1 2000000 3000000" << endl;
    cerr << " chr1 3000000 4000000" << endl;
    cerr << " chr1 4000000 5000000" << endl;
    cerr << " ..." << endl;
    cerr << endl;

    cerr << " # Divide the human genome into sliding (=overlapping) windows of 1MB, with 500KB overlap:" << endl;
    cerr << " $ " << PROGRAM_NAME << " -g hg19.txt -w 1000000 -s 500000" << endl;
    cerr << " chr1 0 1000000" << endl;
    cerr << " chr1 500000 1500000" << endl;
    cerr << " chr1 1000000 2000000" << endl;
    cerr << " chr1 1500000 2500000" << endl;
    cerr << " chr1 2000000 3000000" << endl;
    cerr << " ..." << endl;
    cerr << endl;

    cerr << " # Divide each chromosome in human genome to 1000 windows of equal size:" << endl;
    cerr << " $ " << PROGRAM_NAME << " -g hg19.txt -n 1000" << endl;
    cerr << " chr1 0 249251" << endl;
    cerr << " chr1 249251 498502" << endl;
    cerr << " chr1 498502 747753" << endl;
    cerr << " chr1 747753 997004" << endl;
    cerr << " chr1 997004 1246255" << endl;
    cerr << " ..." << endl;
    cerr << endl;

    cerr << " # Divide each interval in the given BED file into 10 equal-sized windows:" << endl;
    cerr << " $ cat input.bed" << endl;
    cerr << " chr5 60000 70000" << endl;
    cerr << " chr5 73000 90000" << endl;
    cerr << " chr5 100000 101000" << endl;
    cerr << " $ " << PROGRAM_NAME << " -b input.bed -n 10" << endl;
    cerr << " chr5 60000 61000" << endl;
    cerr << " chr5 61000 62000" << endl;
    cerr << " chr5 62000 63000" << endl;
    cerr << " chr5 63000 64000" << endl;
    cerr << " chr5 64000 65000" << endl;
    cerr << " ..." << endl;
    cerr << endl;

    cerr << " # Add a name column, based on the window number: "<< endl;
    cerr << " $ cat input.bed" << endl;
    cerr << " chr5  60000  70000 AAA" << endl;
    cerr << " chr5  73000  90000 BBB" << endl;
    cerr << " chr5 100000 101000 CCC" << endl;
    cerr << " $ " << PROGRAM_NAME << " -b input.bed -n 3 -i winnum" << endl;
    cerr << " chr5        60000   63334   1" << endl;
    cerr << " chr5        63334   66668   2" << endl;
    cerr << " chr5        66668   70000   3" << endl;
    cerr << " chr5        73000   78667   1" << endl;
    cerr << " chr5        78667   84334   2" << endl;
    cerr << " chr5        84334   90000   3" << endl;
    cerr << " chr5        100000  100334  1" << endl;
    cerr << " chr5        100334  100668  2" << endl;
    cerr << " chr5        100668  101000  3" << endl;
    cerr << " ..." << endl;
    cerr << endl;

    cerr << " # Reverse window numbers: "<< endl;
    cerr << " $ cat input.bed" << endl;
    cerr << " chr5  60000  70000 AAA" << endl;
    cerr << " chr5  73000  90000 BBB" << endl;
    cerr << " chr5 100000 101000 CCC" << endl;
    cerr << " $ " << PROGRAM_NAME << " -b input.bed -n 3 -i winnum -reverse" << endl;
    cerr << " chr5        60000   63334   3" << endl;
    cerr << " chr5        63334   66668   2" << endl;
    cerr << " chr5        66668   70000   1" << endl;
    cerr << " chr5        73000   78667   3" << endl;
    cerr << " chr5        78667   84334   2" << endl;
    cerr << " chr5        84334   90000   1" << endl;
    cerr << " chr5        100000  100334  3" << endl;
    cerr << " chr5        100334  100668  2" << endl;
    cerr << " chr5        100668  101000  1" << endl;
    cerr << " ..." << endl;
    cerr << endl;


    cerr << " # Add a name column, based on the source ID + window number: "<< endl;
    cerr << " $ cat input.bed" << endl;
    cerr << " chr5  60000  70000 AAA" << endl;
    cerr << " chr5  73000  90000 BBB" << endl;
    cerr << " chr5 100000 101000 CCC" << endl;
    cerr << " $ " << PROGRAM_NAME << " -b input.bed -n 3 -i srcwinnum" << endl;
    cerr << " chr5        60000   63334   AAA_1" << endl;
    cerr << " chr5        63334   66668   AAA_2" << endl;
    cerr << " chr5        66668   70000   AAA_3" << endl;
    cerr << " chr5        73000   78667   BBB_1" << endl;
    cerr << " chr5        78667   84334   BBB_2" << endl;
    cerr << " chr5        84334   90000   BBB_3" << endl;
    cerr << " chr5        100000  100334  CCC_1" << endl;
    cerr << " chr5        100334  100668  CCC_2" << endl;
    cerr << " chr5        100668  101000  CCC_3" << endl;
    cerr << " ..." << endl;
    cerr << endl;



    cerr << endl;


    exit(1);

}

