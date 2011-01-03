/*****************************************************************************
  bed12ToBed6.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "bedFile.h"
#include "version.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bed12ToBed6"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void ShowHelp(void);
void DetermineBedInput(BedFile *bed);
void ProcessBed(istream &bedInput, BedFile *bed);



int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile       = "stdin";
    bool haveBed         = true;

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
                bedFile = argv[i + 1];
                i++;
            }
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have an input files
    if (!haveBed ) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i (BED) file. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedFile *bed       = new BedFile(bedFile);
        DetermineBedInput(bed);
    }
    else {
        ShowHelp();
    }
}


void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Splits BED12 features into discrete BED6 features." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed12>" << endl << endl;


    // end the program here
    exit(1);
}


void DetermineBedInput(BedFile *bed) {

    // dealing with a proper file
    if (bed->bedFile != "stdin") {

        ifstream bedStream(bed->bedFile.c_str(), ios::in);
        if ( !bedStream ) {
            cerr << "Error: The requested bed file (" << bed->bedFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        ProcessBed(bedStream, bed);
    }
    // reading from stdin
    else {
        ProcessBed(cin, bed);
    }
}


void ProcessBed(istream &bedInput, BedFile *bed) {

    // process each BED entry and convert to BAM
    BED bedEntry, nullBed;
    int lineNum = 0;
    BedLineStatus bedStatus;
    // open the BED file for reading.
    bed->Open();
    while ((bedStatus = bed->GetNextBed(bedEntry, lineNum)) != BED_INVALID) {
        if (bedStatus == BED_VALID) {

            bedVector bedBlocks;  // vec to store the discrete BED "blocks" from a
            splitBedIntoBlocks(bedEntry, lineNum, bedBlocks);

            vector<BED>::const_iterator bedItr  = bedBlocks.begin();
            vector<BED>::const_iterator bedEnd  = bedBlocks.end();
            for (; bedItr != bedEnd; ++bedItr) {
                printf ("%s\t%d\t%d\t%s\t%s\t%s\n", bedItr->chrom.c_str(), bedItr->start, bedItr->end, bedItr->name.c_str(),
                                                    bedItr->score.c_str(), bedItr->strand.c_str());
            }
            bedEntry = nullBed;
        }
    }
    // close up
    bed->Close();
}
