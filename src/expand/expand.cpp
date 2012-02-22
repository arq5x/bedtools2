/*****************************************************************************
expand.cpp

(c) 2009, 2010, 2011 - Aaron Quinlan
Center for Public Health Genomics
University of Virginia
aaronquinlan@gmail.com

Licenced under the MIT license.
******************************************************************************/
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <exception>
#include <stdexcept> // out_of_range exception

#include "version.h"
#include "lineFileUtilities.h"
#include "tabFile.h"
#include "VectorOps.h"
using namespace std;


// define our program name
#define PROGRAM_NAME "bedtools expand"
// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) ((strncmp(argv[i], param, min(actualLen, paramLen))== 0) && \
                                                              (actualLen == paramLen))
#define LOOKS_LIKE_A_PARAM(string) (strlen(string)>0 && string[0]=='-')

// function declarations
void expand_help(void);
void Expand(const string &inFile, 
            const vector<int> &expColumns);

int expand_main(int argc, char* argv[]) {

    // input files
    string inFile             = "stdin";
    string groupColumnsString = "1,2,3";
    string expColumnString;

    // our configuration variables
    bool showHelp          = false;
    bool haveExpColumns     = false;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) expand_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                inFile     = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-c", 2, parameterLength)) {
            if ((i+1) >= argc || LOOKS_LIKE_A_PARAM(argv[i+1])) {
                cerr << endl << "*****ERROR: -opCols parameter requires a value." << endl << endl;
                expand_help();
                break;
            }
            else {
                haveExpColumns       = true;
                expColumnString     = argv[i + 1];
                i++;
            }
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    if (!haveExpColumns) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -opCols." << endl << "*****" << endl;
        showHelp = true;
    }
    


    if (!showHelp) {
        vector<int> expColumns;
        Tokenize(expColumnString, expColumns, ',');

        // sanity check the exp columns
        for(size_t i = 0; i < expColumns.size(); ++i) {
            int expCol = expColumns[i];
            if (expCol < 1) {
                cerr << endl << "*****" << endl << "*****ERROR: expansion columns must be >=1. " << endl << "*****" << endl;
                expand_help();
            }
        }
        Expand(inFile, expColumns);
    }
    else {
        expand_help();
    }
    return 0;
}

void expand_help(void) {

    cerr << "\nTool:    bedtools expand " << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Replicate lines in a file based on columns of comma-separated values." << endl << endl;

    cerr << "Usage:\t " << PROGRAM_NAME << " -c [COLS] " << endl;

    cerr << "Options: " << endl;
    cerr << "\t-i\t"        << "Input file. Assumes \"stdin\" if omitted." << endl << endl;

    cerr << "\t-c \t"    << "Specify the column (1-based) that should be summarized." << endl;
    cerr                 << "\t\t- Required." << endl;

    cerr << "Examples: " << endl;
    cerr << "  $ cat test.txt" << endl;
    cerr << "  chr1	10	20	1,2,3	10,20,30" << endl;
    cerr << "  chr1	40	50	4,5,6	40,50,60" << endl << endl;
               
    cerr << "  $ bedtools expand test.txt -c 5" << endl;
    cerr << "  chr1	10	20	1,2,3	10" << endl;
    cerr << "  chr1	10	20	1,2,3	20" << endl;
    cerr << "  chr1	10	20	1,2,3	30" << endl;
    cerr << "  chr1	40	50	4,5,6	40" << endl;
    cerr << "  chr1	40	50	4,5,6	50" << endl;
    cerr << "  chr1	40	50	4,5,6	60" << endl << endl;

    cerr << "  $ bedtools expand test.txt -c 4,5" << endl;
    cerr << "  chr1	10	20	1	10" << endl;
    cerr << "  chr1	10	20	2	20" << endl;
    cerr << "  chr1	10	20	3	30" << endl;
    cerr << "  chr1	40	50	4	40" << endl;
    cerr << "  chr1	40	50	5	50" << endl;
    cerr << "  chr1	40	50	6	60" << endl;

    // end the program here
    exit(1);

}


void Expand (const string &inFile,
    const vector<int> &expColumns) 
{

    // current line number
    int lineNum = 0;
    // string representing current line
    string inLine;

    // vector of strings holding the tokenized current line
    vector<string>  inFields;
    inFields.reserve(20);

    // build a map of the columns to be expanded
    // to allow quic lookups to test if a column is
    // "normal" or whether it is one of the columns 
    // that is being expaded
    map<int, bool> expColMap;
    for (size_t c = 0; c < expColumns.size(); c++)
        expColMap[expColumns[c]] = true;
    
    // open a new tab file, loop through it line by line
    // and expand each line into multiple lines according to the
    // columns the user has requested.
    // 
    TabLineStatus tabLineStatus;
    TabFile *_tab = new TabFile(inFile);
    _tab->Open();
    while ((tabLineStatus = _tab->GetNextTabLine(inFields, lineNum)) != TAB_INVALID) {
        lineNum++;
        if (tabLineStatus == TAB_VALID) {
            
            // a list containing the expanded values (inner) for each column (outer)
            vector< vector<string> >  expandedCols;
            
            // expand each requested column into a vector
            int prev_size = -1;
            for (size_t c = 0; c < expColumns.size(); c++)
            {
                vector<string> expansion;
                if ((expColumns[c]-1) >= (int) inFields.size()) {
                    cerr << endl
                    << "*****" << endl
                    << "*****ERROR: Requested column number exceeds number of columns." << endl
                    << "*****       This was violated at line: " << lineNum << endl
                    << "*****" << endl;
                    exit(1);
                }
               
                // expand the requested column into a vector
                Tokenize(inFields[expColumns[c]-1], expansion, ',');
               
                if ((int) expansion.size() != prev_size && prev_size >= 0) {
                    cerr << endl
                    << "*****" << endl
                    << "*****ERROR: Each expanded column must have the same number of elements." << endl
                    << "*****       This was violated at line: " << lineNum << endl
                    << "*****" << endl;
                    exit(1);
                }
                else {
                    expandedCols.push_back(expansion);
                }
                prev_size = expansion.size();
            }

            // now replicate/expand the original line based on the 
            // values in the requested columns
            size_t totalCols  = inFields.size();
            for (size_t n = 0; n < expandedCols[0].size(); n++) 
            {
                int numExpColsSeen = 0;
                for (size_t c = 0; c < totalCols; c++)
                {
                    // normal column, print as-is
                    if (!expColMap[c+1]) {
                        printf("%s", inFields[c].c_str());
                    }
                    // expanded column, grab relevant value from expanded vector.
                    else {
                        // expandedCols[i][j]
                        // i == column, j = row
                        printf("%s", expandedCols[numExpColsSeen][n].c_str());
                        numExpColsSeen++;
                    }
                    // add a tab if not the very last value
                    if (c < totalCols - 1)
                        printf("\t");
                }
                printf("\n");
            }
        }
        inFields.clear();
    }
    _tab->Close();
}

