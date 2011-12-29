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

#include "genomeFile.h"
#include "multiIntersectBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools multiinter"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

//STLized version of basename()
// (because POSIX basename() modifies the input string pointer)
// Additionally: removes any extension the basename might have.
std::string stl_basename(const std::string& path);

// function declarations
void multiintersect_help(void);
void multiintersect_examples(void);

int multiintersect_main(int argc, char* argv[])
{
    bool haveFiles         = false;
    bool haveTitles        = false;
    bool haveGenome        = false;
    bool haveFiller        = true;
    bool printHeader       = false;
    bool printEmptyRegions = false;
    bool cluster           = false;
    bool showHelp          = false;
    string genomeFile;
    string basePath;
    string noCoverageValue = "0";
    vector<string> inputFiles;
    vector<string> inputTitles;

    //Parse command line options
    if(argc <= 1)
        multiintersect_help();

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp == true) {
        multiintersect_help();
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
        else if(PARAMETER_CHECK("-cluster", 8, parameterLength)) {
            cluster = true;
        }
        else if(PARAMETER_CHECK("-examples", 9, parameterLength)) {
            multiintersect_help();
            multiintersect_examples();
            exit(1);
        }
    }

    //Sanity checks
    if (inputFiles.empty() == true) {
        cerr << "Error: missing file names (-i) to combine." << endl;
        exit(1);
    }
    if (inputFiles.size() == 1) {
        cerr << "Error: Only a single file was specified. Nothing to combine, exiting." << endl;
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

    MultiIntersectBed mbi(cout, inputFiles, inputTitles, printEmptyRegions, genomeFile, noCoverageValue);
    if (printHeader)
        mbi.PrintHeader();
    if (!cluster)
        mbi.MultiIntersect();
    else
        mbi.Cluster();
    
    return 0;
}

void multiintersect_help(void) {

    cerr << "\nTool:    bedtools multiinter (aka multiIntersectBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Identifies common intervals among multiple" << endl;
    cerr << "\t BED/GFF/VCF files." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i FILE1 FILE2 .. FILEn" << endl;
    cerr << "\t Requires that each interval file is sorted by chrom/start. " << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-cluster\t"      << "Invoke Ryan Layers's clustering algorithm." << endl << endl;

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



void multiintersect_examples()
{
    cerr << "Example usage:\n\n"  \
"== Input files: ==\n" \
"\n" \
" $ cat a.bed\n" \
" chr1  6   12\n" \
" chr1  10  20\n" \
" chr1  22  27\n" \
" chr1  24  30\n" \
"\n" \
" $ cat b.bed\n" \
" chr1  12  32\n" \
" chr1  14  30\n" \
"\n" \
" $ cat c.bed\n" \
" chr1  8   15\n" \
" chr1  10  14\n" \
" chr1  32  34\n" \
"\n" \
" $ cat sizes.txt\n" \
" chr1  5000\n" \
"\n" \
"== Multi-intersect the files: ==\n" \
"\n" \
" $ multiIntersectBed -i a.bed b.bed c.bed\n" \
"chr1	6	8	1	1	1	0	0\n" \
"chr1	8	12	2	1,3	1	0	1\n" \
"chr1	12	15	3	1,2,3	1	1	1\n" \
"chr1	15	20	2	1,2	1	1	0\n" \
"chr1	20	22	1	2	0	1	0\n" \
"chr1	22	30	2	1,2	1	1	0\n" \
"chr1	30	32	1	2	0	1	0\n" \
"chr1	32	34	1	3	0	0	1\n" \
"\n" \
"== Multi-intersect the files, with a header line (titles are the file names): ==\n" \
"\n" \
" $ multiIntersectBed -header -i a.bed b.bed c.bed\n" \
" chrom	start	end	num	list	a.bed	b.bed	c.bed\n" \
" chr1	6	8	1	1	1	0	0\n" \
" chr1	8	12	2	1,3	1	0	1\n" \
" chr1	12	15	3	1,2,3	1	1	1\n" \
" chr1	15	20	2	1,2	1	1	0\n" \
" chr1	20	22	1	2	0	1	0\n" \
" chr1	22	30	2	1,2	1	1	0\n" \
" chr1	30	32	1	2	0	1	0\n" \
" chr1	32	34	1	3	0	0	1\n" \
"\n" \
"== Multi-intersect the files, with a header line and custom names: ==\n" \
"\n" \
" $ multiIntersectBed -header -i a.bed b.bed c.bed -names A B C\n" \
" chrom	start	end	num	list	A	B	C\n" \
" chr1	6	8	1	A	1	0	0\n" \
" chr1	8	12	2	A,C	1	0	1\n" \
" chr1	12	15	3	A,B,C	1	1	1\n" \
" chr1	15	20	2	A,B	1	1	0\n" \
" chr1	20	22	1	B	0	1	0\n" \
" chr1	22	30	2	A,B	1	1	0\n" \
" chr1	30	32	1	B	0	1	0\n" \
" chr1	32	34	1	C	0	0	1\n" \
"\n" \
"== Multi-intersect the files, showing empty regions (note, requires -g): ==\n" \
"\n" \
" $ multiIntersectBed -header -i a.bed b.bed c.bed -names A B C -empty -g sizes.txt\n" \
" chrom	start	end	num	list	A	B	C\n" \
" chr1	0	6	0	none	0	0	0\n" \
" chr1	6	8	1	A	1	0	0\n" \
" chr1	8	12	2	A,C	1	0	1\n" \
" chr1	12	15	3	A,B,C	1	1	1\n" \
" chr1	15	20	2	A,B	1	1	0\n" \
" chr1	20	22	1	B	0	1	0\n" \
" chr1	22	30	2	A,B	1	1	0\n" \
" chr1	30	32	1	B	0	1	0\n" \
" chr1	32	34	1	C	0	0	1\n" \
" chr1	34	5000	0	none	0	0	0\n" \
"\n"
;
}

std::string stl_basename(const std::string& path)
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
