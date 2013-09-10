/*****************************************************************************
  unionBedGraphsMain.cpp

  (c) 2010 - Assaf Gordon, CSHL
           - Aaron Quinlan, UVA
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <climits>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <getopt.h>
#include <libgen.h> //for basename()

#include "GenomeFile.h"
#include "unionBedGraphs.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools unionbedg"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

//STLized version of basename()
// (because POSIX basename() modifies the input string pointer)
// Additionally: removes any extension the basename might have.
std::string ubg_stl_basename(const std::string& path);

// function declarations
void unionbedgraphs_help(void);
void unionbedgraphs_showexamples(void);


int unionbedgraphs_main(int argc, char* argv[])
{
    bool haveFiles         = false;
    bool haveTitles        = false;
    bool haveGenome        = false;
    bool haveFiller        = true;
    bool printHeader       = false;
    bool printEmptyRegions = false;
    bool showHelp          = false;
    string genomeFile;
    string basePath;
    string noCoverageValue = "0";
    vector<string> inputFiles;
    vector<string> inputTitles;

    //Parse command line options
    if(argc <= 1)
        unionbedgraphs_help();

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp == true) {
        unionbedgraphs_help();
        exit(1);
    }

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveFiles = true;
                i = i+1;
                string file = argv[i];
                while (file[0] != '-' && i < argc) {
                    inputFiles.push_back(file);
                    i++;
                    if (i < argc)
                        file = argv[i];
                }
                i--;
            }
        }
        else if(PARAMETER_CHECK("-names", 6, parameterLength)) {
            if ((i+1) < argc) {
                haveTitles = true;
                i = i+1;
                string title = argv[i];
                while (title[0] != '-' && i < argc) {
                    inputTitles.push_back(title);
                    i++;
                    if (i < argc)
                        title = argv[i];
                }
                i--;
            }
        }
        else if(PARAMETER_CHECK("-g", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveGenome = true;
                genomeFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-filler", 7, parameterLength)) {
            if ((i+1) < argc) {
                haveFiller      = true;
                noCoverageValue = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-header", 7, parameterLength)) {
            printHeader = true;
        }
        else if(PARAMETER_CHECK("-empty", 6, parameterLength)) {
            printEmptyRegions = true;
        }
        else if(PARAMETER_CHECK("-examples", 9, parameterLength)) {
            unionbedgraphs_help();
            unionbedgraphs_showexamples();
            exit(1);
        }
    }

    //Sanity checks
    if (inputFiles.empty() == true) {
        cerr << "Error: missing BedGraph file names (-i) to combine." << endl;
        exit(1);
    }
    if (inputFiles.size() == 1) {
        cerr << "Error: Only a single BedGraph file was specified. Nothing to combine, exiting." << endl;
        exit(1);
    }
    if (printEmptyRegions && (genomeFile.empty() == true)) {
        cerr << "Error: when using -empty, the genome sizes file (-g) must be specified using '-g FILE'." << endl;
        exit(1);
    }
    if ((haveTitles == true) && (inputFiles.size() != inputTitles.size())) {
        cerr << "Error: The number of file titles (-names) does not match the number of files (-i)." << endl;
        exit(1);
    }

    UnionBedGraphs ubg(cout, inputFiles, inputTitles, printEmptyRegions, genomeFile, noCoverageValue);
    if (printHeader)
        ubg.PrintHeader();
    ubg.Union();
    
    return 0;
}

void unionbedgraphs_help(void) {

    cerr << "\nTool:    bedtools unionbedg (aka unionBedGraphs)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Combines multiple BedGraph files into a single file," << endl;
    cerr << "\t allowing coverage comparisons between them." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i FILE1 FILE2 .. FILEn" << endl;
    cerr << "\t Assumes that each BedGraph file is sorted by chrom/start " << endl;
    cerr << "\t and that the intervals in each are non-overlapping." << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-header\t\t"     << "Print a header line." << endl;
    cerr                        << "\t\t\t(chrom/start/end + names of each file)." << endl << endl;

    cerr << "\t-names\t\t"      << "A list of names (one/file) to describe each file in -i." << endl;
    cerr                        << "\t\t\tThese names will be printed in the header line." << endl << endl;

    cerr << "\t-g\t\t"          << "Use genome file to calculate empty regions." << endl;
    cerr                        << "\t\t\t- STRING." << endl << endl;

    cerr << "\t-empty\t\t"      << "Report empty regions (i.e., start/end intervals w/o" << endl;
    cerr                        << "\t\t\tvalues in all files)." << endl;
    cerr                        << "\t\t\t- Requires the '-g FILE' parameter.\n" << endl;

    cerr << "\t-filler TEXT\t"  << "Use TEXT when representing intervals having no value." << endl;
    cerr                        << "\t\t\t- Default is '0', but you can use 'N/A' or any text." << endl << endl;

    cerr << "\t-examples\t"     << "Show detailed usage examples." << endl << endl;
}



void unionbedgraphs_showexamples()
{
    cerr << "Example usage:\n\n"  \
"== Input files: ==\n" \
"\n" \
" $ cat 1.bg\n" \
" chr1  1000    1500    10\n" \
" chr1  2000    2100    20\n" \
"\n" \
" $ cat 2.bg\n" \
" chr1  900 1600    60\n" \
" chr1  1700    2050    50\n" \
"\n" \
" $ cat 3.bg\n" \
" chr1  1980    2070    80\n" \
" chr1  2090    2100    20\n" \
"\n" \
" $ cat sizes.txt\n" \
" chr1  5000\n" \
"\n" \
"== Union/combine the files: ==\n" \
"\n" \
" $ unionBedGraphs -i 1.bg 2.bg 3.bg\n" \
" chr1  900 1000    0   60  0\n" \
" chr1  1000    1500    10  60  0\n" \
" chr1  1500    1600    0   60  0\n" \
" chr1  1700    1980    0   50  0\n" \
" chr1  1980    2000    0   50  80\n" \
" chr1  2000    2050    20  50  80\n" \
" chr1  2050    2070    20  0   80\n" \
" chr1  2070    2090    20  0   0\n" \
" chr1  2090    2100    20  0   20\n" \
"\n" \
"== Union/combine the files, with a header line (titles are the file names): ==\n" \
"\n" \
" $ unionBedGraphs -header -i 1.bg 2.bg 3.bg\n" \
" chrom start   end 1   2   3\n" \
" chr1  900 1000    0   60  0\n" \
" chr1  1000    1500    10  60  0\n" \
" chr1  1500    1600    0   60  0\n" \
" chr1  1700    1980    0   50  0\n" \
" chr1  1980    2000    0   50  80\n" \
" chr1  2000    2050    20  50  80\n" \
" chr1  2050    2070    20  0   80\n" \
" chr1  2070    2090    20  0   0\n" \
" chr1  2090    2100    20  0   20\n" \
"\n" \
"== Union/combine the files, with a header line and custom names: ==\n" \
"\n" \
" $ unionBedGraphs -header -i 1.bg 2.bg 3.bg -names WT-1 WT-2 KO-1\n" \
" chrom start   end WT-1    WT-2    KO-1\n" \
" chr1  900 1000    0   60  0\n" \
" chr1  1000    1500    10  60  0\n" \
" chr1  1500    1600    0   60  0\n" \
" chr1  1700    1980    0   50  0\n" \
" chr1  1980    2000    0   50  80\n" \
" chr1  2000    2050    20  50  80\n" \
" chr1  2050    2070    20  0   80\n" \
" chr1  2070    2090    20  0   0\n" \
" chr1  2090    2100    20  0   20\n" \
"\n" \
"== Union/combine, showing empty regions (note, requires -g): ==\n" \
"\n" \
" $ unionBedGraphs -header -empty -g sizes.TXT -i 1.bg 2.bg 3.bg\n" \
" chrom start   end 1   2   3\n" \
" chr1  0   900 0   0   0\n" \
" chr1  900 1000    0   60  0\n" \
" chr1  1000    1500    10  60  0\n" \
" chr1  1500    1600    0   60  0\n" \
" chr1  1600    1700    0   0   0\n" \
" chr1  1700    1980    0   50  0\n" \
" chr1  1980    2000    0   50  80\n" \
" chr1  2000    2050    20  50  80\n" \
" chr1  2050    2070    20  0   80\n" \
" chr1  2070    2090    20  0   0\n" \
" chr1  2090    2100    20  0   20\n" \
" chr1  2100    5000    0   0   0\n" \
"\n" \
;
}

std::string ubg_stl_basename(const std::string& path)
{
    string result;

    char* path_dup = strdup(path.c_str());
    char* basename_part = basename(path_dup);
    result = basename_part;
    free(path_dup);

    size_t pos = result.find_last_of('.');
    if (pos != string::npos )
        result = result.substr(0,pos);

    return result;
}

