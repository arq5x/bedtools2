/*****************************************************************************
  sortBedMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "sortBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools sort"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void sort_help(void);

int sort_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile;
    string faidxFile; 
    bool haveBed    = false;
    int sortChoices = 0;

    bool sortBySizeAsc            = false;
    bool sortBySizeDesc           = false;
    bool sortByChromThenSizeAsc   = false;
    bool sortByChromThenSizeDesc  = false;
    bool sortByChromThenScoreAsc  = false;
    bool sortByChromThenScoreDesc = false;
    bool sortByFaidx              = false;
    bool printHeader        = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) sort_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedFile = argv[i + 1];
                i++;
                haveBed = true;
            }
        }
        else if(PARAMETER_CHECK("-sizeA", 6, parameterLength)) {
            sortBySizeAsc = true;
            sortChoices++;
        }
        else if(PARAMETER_CHECK("-sizeD", 6, parameterLength)) {
            sortBySizeDesc = true;
            sortChoices++;
        }
        else if(PARAMETER_CHECK("-chrThenSizeA", 13, parameterLength)) {
            sortByChromThenSizeAsc = true;
            sortChoices++;
        }
        else if(PARAMETER_CHECK("-chrThenSizeD", 13, parameterLength)) {
            sortByChromThenSizeDesc = true;
            sortChoices++;
        }
        else if(PARAMETER_CHECK("-chrThenScoreA", 14, parameterLength)) {
            sortByChromThenScoreAsc = true;
            sortChoices++;
        }
        else if(PARAMETER_CHECK("-chrThenScoreD", 14, parameterLength)) {
            sortByChromThenScoreDesc = true;
            sortChoices++;
        }
        else if(PARAMETER_CHECK("-faidx", 6, parameterLength)) {
             if ((i+1) < argc) {
                faidxFile = argv[i + 1];
                i++;
           		  }
           	sortByFaidx = true;
            sortChoices++;
        }
        else if(PARAMETER_CHECK("-g", 2, parameterLength)) {
             if ((i+1) < argc) {
                faidxFile = argv[i + 1];
                i++;
                  }
            sortByFaidx = true;
            sortChoices++;
        }
        else if(PARAMETER_CHECK("-header", 7, parameterLength)) {
            printHeader = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have input 
    if (!haveBed) {
        if (!isatty(STDIN_FILENO))
        {
            bedFile = "stdin";
        }
        else 
        {
            cerr << endl << "*****" << endl << "*****ERROR: Need -i BED file. " << endl << "*****" << endl;
            showHelp = true;
        }
    }
    if (sortChoices > 1) {
        cerr << endl << "*****" << endl << "*****ERROR: Sorting options are mutually exclusive.  Please choose just one. " << endl << "*****" << endl;
        showHelp = true;
    }


    if (!showHelp) {
        BedSort *bm = new BedSort(bedFile, printHeader,faidxFile);

        if (sortBySizeAsc) {
            bm->SortBedBySizeAsc();
        }
        else if (sortBySizeDesc) {
            bm->SortBedBySizeDesc();
        }
        else if (sortByChromThenSizeAsc) {
            bm->SortBedByChromThenSizeAsc();
        }
        else if (sortByChromThenSizeDesc) {
            bm->SortBedByChromThenSizeDesc();
        }
        else if (sortByChromThenScoreAsc) {
            bm->SortBedByChromThenScoreAsc();
        }
        else if (sortByChromThenScoreDesc) {
            bm->SortBedByChromThenScoreDesc();
        }
        else if (sortByFaidx) {
            bm->SortBedOnFaidx();
        }       
        else {
            bm->SortBed();
        }
        return 0;
    }
    else {
        sort_help();
    }
    return 0;
}

void sort_help(void) {

    cerr << "\nTool:    bedtools sort (aka sortBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Sorts a feature file in various and useful ways." << endl << endl;
    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t" << "-sizeA\t\t\t"    << "Sort by feature size in ascending order." << endl;
    cerr << "\t" << "-sizeD\t\t\t"    << "Sort by feature size in descending order." << endl;
    cerr << "\t" << "-chrThenSizeA\t\t"   << "Sort by chrom (asc), then feature size (asc)." << endl;
    cerr << "\t" << "-chrThenSizeD\t\t"   << "Sort by chrom (asc), then feature size (desc)." << endl;
    cerr << "\t" << "-chrThenScoreA\t\t"  << "Sort by chrom (asc), then score (asc)." << endl;
    cerr << "\t" << "-chrThenScoreD\t\t"  << "Sort by chrom (asc), then score (desc)." << endl;
    cerr << "\t" << "-g (names.txt)\t"    << "Sort according to the chromosomes declared in \"genome.txt\"" << endl;
    cerr << "\t" << "-faidx (names.txt)\t"  << "Sort according to the chromosomes declared in \"names.txt\"" << endl;
    
    cerr << "\t-header\t"       << "Print the header from the A file prior to results." << endl << endl;

    exit(1);

}
