/*****************************************************************************
  overlap.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "lineFileUtilities.h"
#include "bedFile.h"
#include "version.h"
using namespace std;


// define our program name
#define PROGRAM_NAME "getOverlap"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void getoverlap_help(void);
void DetermineInput(string &inFile, short &s1Col, short &e1Col, short &s2Col, short &e2Col);
void ComputeOverlaps(istream &input, short &s1Col, short &e1Col, short &s2Col, short &e2Col);

int getoverlap_main(int argc, char* argv[]) {

    // input files
    string inFile = "stdin";
    string columns;

    // our configuration variables
    bool showHelp = false;
    bool haveInFile  = true;
    bool haveColumns = false;


    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) getoverlap_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                inFile     = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-cols", 5, parameterLength)) {
            haveColumns = true;
            columns     = argv[i + 1];
            i++;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have an input files
    if (!haveInFile ) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i file. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {

        // Split the column string sent by the user into discrete column numbers
        // A comma separated string is expected.
        vector<string> posColumns;
        Tokenize(columns, posColumns, ',');

        if (posColumns.size() != 4) {
            cerr << endl << "*****" << endl << "*****ERROR: Please specify 4, comma-separated position columns. " << endl << "*****" << endl;
            getoverlap_help();
        }
        else {
            short s1, e1, s2, e2;
            s1 = atoi(posColumns[0].c_str());
            e1 = atoi(posColumns[1].c_str());
            s2 = atoi(posColumns[2].c_str());
            e2 = atoi(posColumns[3].c_str());

            DetermineInput(inFile, s1, e1, s2, e2);
        }
    }
    else {
        getoverlap_help();
    }
    return 0;
}

void getoverlap_help(void) {

    cerr << "\nTool:    bedtools overlap (aka getOverlap)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Computes the amount of overlap (positive values)" << endl;
    cerr << "\t or distance (negative values) between genome features" << endl;
    cerr << "\t and reports the result at the end of the same line." << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-i\t"        << "Input file. Use \"stdin\" for pipes." << endl << endl;

    cerr << "\t-cols\t"     << "Specify the columns (1-based) for the starts and ends of the" << endl;
    cerr                    << "\t\tfeatures for which you'd like to compute the overlap/distance." << endl;
    cerr                    << "\t\tThe columns must be listed in the following order: " << endl << endl;
    cerr                    << "\t\tstart1,end1,start2,end2" << endl << endl;

    cerr << "Example: " << endl;
    cerr << "\t$ windowBed -a A.bed -b B.bed -w 10" << endl;
    cerr << "\tchr1 10  20  A   chr1    15  25  B" << endl;
    cerr << "\tchr1 10  20  C   chr1    25  35  D" << endl << endl;
    cerr << "\t$ windowBed -a A.bed -b B.bed -w 10 | overlap -i stdin -cols 2,3,6,7" << endl;
    cerr << "\tchr1 10  20  A   chr1    15  25  B   5" << endl;
    cerr << "\tchr1 10  20  C   chr1    25  35  D   -5" << endl;

    // end the program here
    exit(1);

}


void DetermineInput(string &inFile, short &s1Col, short &e1Col, short &s2Col, short &e2Col) {


    if (inFile != "stdin") {   // process a file

        ifstream in(inFile.c_str(), ios::in);
        if ( !in ) {
            cerr << "Error: The requested input file (" << inFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        ComputeOverlaps(in, s1Col, e1Col, s2Col, e2Col);
    }
    else ComputeOverlaps(cin, s1Col, e1Col, s2Col, e2Col);
}


void ComputeOverlaps(istream &input, short &s1Col, short &e1Col, short &s2Col, short &e2Col) {

    int lineNum = 0;
    string inLine;
    vector<string> inFields;

    int overlap;

    char *s1End, *e1End, *s2End, *e2End;
    long s1, e1, s2, e2;

    while (getline(input, inLine)) {
        lineNum++;
        Tokenize(inLine, inFields);

        if (inFields.size() > 1) {

            // test if columns  2 and 3 are integers.  If so, assume BED.
            s1 = strtol(inFields[s1Col-1].c_str(), &s1End, 10);
            e1 = strtol(inFields[e1Col-1].c_str(), &e1End, 10);
            s2 = strtol(inFields[s2Col-1].c_str(), &s2End, 10);
            e2 = strtol(inFields[e2Col-1].c_str(), &e2End, 10);

            // strtol will set pointers to the start of the string if non-integral, base 10
            // if they all check out, we have valid numeric columns.  Otherwise, complain.
            if (s1End != inFields[s1Col-1].c_str() &&
                e1End != inFields[e1Col-1].c_str() &&
                s2End != inFields[s2Col-1].c_str() &&
                e2End != inFields[e2Col-1].c_str()) {

                overlap = overlaps(s1, e1, s2, e2);
                printf("%s\t%d\n", inLine.c_str(), overlap);
            }
            else {
                cerr << "One of your columns appears to be non-numeric at line " << lineNum << ". Exiting..." << endl << endl;
                exit(1);
            }
        }
        inFields.clear();
    }
}
