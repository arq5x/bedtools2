/*****************************************************************************
  bedToIgv.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "bedFile.h"
#include "GenomeFile.h"
#include "version.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools igv"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void bedtoigv_help(void);

void DetermineBedInput(BedFile *bed, string path, string sortType, string session,
                        bool collapse, bool useNames, string imageType, int slop);
void ProcessBed(istream &bedInput, BedFile *bed, string path, string sortType, string session,
                        bool collapse, bool useNames, string imageType, int slop);


int bedtoigv_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile   = "stdin";
    string imagePath = "./";
    string sortType  = "none";
    string session   = "none";
    int slop         = 0;
    string imageType = "png";

    bool haveBed         = true;
    bool collapse        = false;
    bool useNames        = false;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) bedtoigv_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-path", 5, parameterLength)) {
            if ((i+1) < argc) {
                imagePath = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-sort", 5, parameterLength)) {
            if ((i+1) < argc) {
                sortType = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-sess", 5, parameterLength)) {
            if ((i+1) < argc) {
                session = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-clps", 5, parameterLength)) {
            collapse = true;
        }
        else if(PARAMETER_CHECK("-name", 5, parameterLength)) {
            useNames = true;
        }
        else if(PARAMETER_CHECK("-slop", 5, parameterLength)) {
            if ((i+1) < argc) {
                slop = atoi(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-img", 4, parameterLength)) {
            if ((i+1) < argc) {
                imageType = argv[i + 1];
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
    if (sortType != "none") {
        if ((sortType != "base")    && (sortType != "position") && (sortType != "strand") &&
            (sortType != "quality") && (sortType != "sample")   && (sortType != "readGroup")) {
                cerr << endl << "*****" << endl << "*****ERROR: Invalid sort option. " << endl << "*****" << endl;
                showHelp = true;
            }
    }
    if (slop < 0) {
        cerr << endl << "*****" << endl << "*****ERROR: Slop must be >= 0. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedFile *bed       = new BedFile(bedFile);
        DetermineBedInput(bed, imagePath, sortType, session, collapse, useNames, imageType, slop);
    }
    else {
        bedtoigv_help();
    }
    return 0;
}


void bedtoigv_help(void) {

    cerr << "\nTool:    bedtools igv (aka bedToIgv)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Creates a batch script to create IGV images " << endl;
    cerr << "         at each interval defined in a BED/GFF/VCF file." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-path\t"     << "The full path to which the IGV snapshots should be written." << endl;
    cerr                    << "\t\t(STRING) Default: ./" << endl << endl;

    cerr << "\t-sess\t"     << "The full path to an existing IGV session file to be " << endl;
    cerr                    << "\t\tloaded prior to taking snapshots." << endl << endl;
    cerr                    << "\t\t(STRING) Default is for no session to be loaded." << endl << endl;

    cerr << "\t-sort\t"     << "The type of BAM sorting you would like to apply to each image. " << endl;
    cerr                    << "\t\tOptions: base, position, strand, quality, sample, and readGroup" << endl;
    cerr                    << "\t\tDefault is to apply no sorting at all." << endl << endl;

    cerr << "\t-clps\t"     << "Collapse the aligned reads prior to taking a snapshot. " << endl;
    cerr                    << "\t\tDefault is to no collapse." << endl << endl;

    cerr << "\t-name\t"     << "Use the \"name\" field (column 4) for each image's filename. " << endl;
    cerr                    << "\t\tDefault is to use the \"chr:start-pos.ext\"." << endl << endl;

    cerr << "\t-slop\t"     << "Number of flanking base pairs on the left & right of the image." << endl;
    cerr                    << "\t\t- (INT) Default = 0." << endl << endl;

    cerr << "\t-img\t"      << "The type of image to be created. " << endl;
    cerr                    << "\t\tOptions: png, eps, svg" << endl;
    cerr                    << "\t\tDefault is png." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  The resulting script is meant to be run from within IGV." << endl;
    cerr << "\t(2)  Unless you use the -sess option, it is assumed that prior to " << endl;
    cerr << "\t\trunning the script, you've loaded the proper genome and tracks." << endl << endl;


    // end the program here
    exit(1);
}


void DetermineBedInput(BedFile *bed, string path, string sortType, string session,
                       bool collapse, bool useNames, string imageType, int slop) {

    // dealing with a proper file
    if (bed->bedFile != "stdin") {

        ifstream bedStream(bed->bedFile.c_str(), ios::in);
        if ( !bedStream ) {
            cerr << "Error: The requested bed file (" << bed->bedFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        ProcessBed(bedStream, bed, path, sortType, session, collapse, useNames, imageType, slop);
    }
    // reading from stdin
    else {
        ProcessBed(cin, bed, path, sortType, session, collapse, useNames, imageType, slop);
    }
}


void ProcessBed(istream &bedInput, BedFile *bed, string path, string sortType, string session,
                bool collapse, bool useNames, string imageType, int slop) {

    // set the image path
    cout << "snapshotDirectory " << path << endl;

    // should we load a session
    if (session != "none")
        cout << "load " << session << endl;


    BED bedEntry;
    bed->Open();
    // process each BED entry and convert to an IGV request
    while (bed->GetNextBed(bedEntry)) {
        if (bed->_status == BED_VALID) {

            string filename = bedEntry.chrom + "_" + ToString(bedEntry.start) + "_" + ToString(bedEntry.end);
            string locus    = bedEntry.chrom + ":" + ToString(bedEntry.start - slop) + "-" + ToString(bedEntry.end + slop);

            if (useNames == true) {
                if (bedEntry.name.empty() == false)
                    filename = filename + "_" + bedEntry.name;
                else {
                    cerr << "Error: You requested that filenames be based upon the name field.  However, it appears to be empty. Exiting!" << endl;
                    exit (1);
                }
            }
            if (slop > 0) {
                filename = filename + "_" + "slop" + ToString(slop);
            }
            // goto
            cout << "goto " << locus << endl;

            // sort
            if (sortType != "none")
                cout << "sort " << sortType << endl;

            // collapse
            if (collapse == true)
                cout << "collapse" << endl;

            // snapshot
            cout << "snapshot " << filename << "." << imageType << endl;

        }
    }
    // close up
    bed->Close();
}
