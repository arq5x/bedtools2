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
#include "BlockedIntervals.h"
#include "bedFile.h"
#include "version.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bedtools bed12tobed6"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void bed12tobed6_help(void);
void DetermineBedInput(BedFile *bed);
void ProcessBed(istream &bedInput, BedFile *bed);


bool addBlockNums = false;

int bed12tobed6_main(int argc, char* argv[]) {

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

    if(showHelp) bed12tobed6_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-n", 2, parameterLength)) {
            addBlockNums = true;
            i++;
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
        bed12tobed6_help();
    }
    return 0;
}


void bed12tobed6_help(void) {

    cerr << "\nTool:    bedtools bed12tobed6 (aka bed12ToBed6)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Splits BED12 features into discrete BED6 features." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed12>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-n\t"        << "Force the score to be the (1-based) block number from the BED12." << endl << endl;


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
    BED bedEntry;
    // open the BED file for reading.
    bed->Open();
    while (bed->GetNextBed(bedEntry)) {
        if (bed->_status == BED_VALID) {

            bedVector bedBlocks;  // vec to store the discrete BED "blocks" from a
            GetBedBlocks(bedEntry, bedBlocks);

            for (int i = 0; i < (int) bedBlocks.size(); ++i) {
                if (addBlockNums == false) {
                    printf ("%s\t%d\t%d\t%s\t%s\t%s\n", bedBlocks[i].chrom.c_str(), bedBlocks[i].start, bedBlocks[i].end, bedBlocks[i].name.c_str(),
                                                        bedBlocks[i].score.c_str(), bedBlocks[i].strand.c_str());
                }
                else {
                    if (bedBlocks[i].strand == "+")
                        printf ("%s\t%d\t%d\t%s\t%d\t%s\n", bedBlocks[i].chrom.c_str(), bedBlocks[i].start, bedBlocks[i].end, bedBlocks[i].name.c_str(),
                                                        i+1, bedBlocks[i].strand.c_str());
                    else 
                        printf ("%s\t%d\t%d\t%s\t%d\t%s\n", bedBlocks[i].chrom.c_str(), bedBlocks[i].start, bedBlocks[i].end, bedBlocks[i].name.c_str(),
                                                        (int) ((bedBlocks.size())-i), bedBlocks[i].strand.c_str());
                }
            }
        }
    }
    // close up
    bed->Close();
}
