/*****************************************************************************
  bedpeIntersect.cpp

  (c) 2025 - Martin Pollard
  Campbell Group
  Genome Research Ltd t/a Wellcome Sanger Institute
  mp15@sanger.ac.uk

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "bedFilePE.h"
#include "GenomeFile.h"
#include "version.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bedpeintersect"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
static void bedpeintersect_help(void);
static void ProcessBedPEs(BedFilePE *bedpe_a, BedFilePE *bedpe_b);
inline bool operator<(BEDPE bedpeEntry_a, BEDPE bedpeEntry_b);
inline bool operator>(BEDPE bedpeEntry_a, BEDPE bedpeEntry_b);
inline bool operator==(BEDPE bedpeEntry_a, BEDPE bedpeEntry_b);

int bedpeintersect_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    std::string bedpeFileA = "stdin";
    std::string bedpeFileB;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) bedpeintersect_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-a", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedpeFileA = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedpeFileB = argv[i + 1];
                i++;
            }
        }
        else {
            std::cerr << std::endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << std::endl << std::endl;
            showHelp = true;
        }
    }

    // make sure we have an input files
    if (bedpeFileB.empty() ) {
        std::cerr << std::endl << "*****" << std::endl << "*****ERROR: Need both -a (BEDPE) and -b (BEDPE) file. " << std::endl << "*****" << std::endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedFilePE *bedpe_a= new BedFilePE(bedpeFileA);
        BedFilePE *bedpe_b= new BedFilePE(bedpeFileB);

        ProcessBedPEs(bedpe_a, bedpe_b);
    }
    else {
        bedpeintersect_help();
    }
    return 0;
}


static void bedpeintersect_help(void) {
    
    std::cerr << "\nTool:    bedtools bedpeintersect (aka bedpeIntersect)" << std::endl;
    std::cerr << "Version: " << VERSION << std::endl;
    std::cerr << "Summary: Forms the intersection of two BEDPE files." << std::endl << std::endl;

    std::cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bedpe> -b <bedpe>" << endl << endl;

    std::cerr << "Options: " << endl;

    // end the program here
    exit(1);
}

static void ProcessBedPEs(BedFilePE *bedpe_a, BedFilePE *bedpe_b) {
    // process each BEDPE entry and compare to B
    BEDPE bedpeEntry_a;
    int lineNum_a = 0;
    BEDPE bedpeEntry_b;
    int lineNum_b = 0;
    BedLineStatus bedpeStatus_a;
    BedLineStatus bedpeStatus_b;

    // open the BEDPE file for reading.
    bedpe_a->Open();
    bedpe_b->Open();
    // Prime the pump
    if ((bedpeStatus_a = bedpe_a->GetNextBedPE(bedpeEntry_a, lineNum_a)) == BED_INVALID) goto ProcessBedPEs_bail;
    bedpeStatus_b = bedpe_b->GetNextBedPE(bedpeEntry_b, lineNum_b);

    // Assume sorted
    do {
        if (bedpeStatus_a == BED_VALID) {
            // Don't read b if it's already EOF
            if (bedpeStatus_b != BED_INVALID) {
                // Read b file until matching coord or past
                do {
                    if (bedpeStatus_b == BED_VALID) {
                        // compare to a are past it?
                        if ( bedpeEntry_a < bedpeEntry_b )
                            break; // if so break out of this loop we need another record from a
                        // If not past it do the comparison
                        if (bedpeEntry_a == bedpeEntry_b) {
                            bedpe_a->reportBedPENewLine(bedpeEntry_a);
                        }
                    }
                } while ((bedpeStatus_b = bedpe_b->GetNextBedPE(bedpeEntry_b, lineNum_b)) != BED_INVALID);
            }
        }
    } while ((bedpeStatus_a = bedpe_a->GetNextBedPE(bedpeEntry_a, lineNum_a)) != BED_INVALID);

ProcessBedPEs_bail:
    //close up
    bedpe_b->Close();
    bedpe_a->Close();
}

inline bool operator<(BEDPE bedpeEntry_a, BEDPE bedpeEntry_b) {
    // TODO: simplify this
    return (bedpeEntry_a.chrom1.compare(bedpeEntry_b.chrom1) <= 0) // Compare chrom1 and chrom1 first
        || ((bedpeEntry_a.chrom1.compare(bedpeEntry_b.chrom1) == 0) && (bedpeEntry_a.start1 < bedpeEntry_b.start1 || (bedpeEntry_a.start1 == bedpeEntry_b.start1 && bedpeEntry_a.end1 < bedpeEntry_b.end1)))
        || ((bedpeEntry_a.chrom1.compare(bedpeEntry_b.chrom1) == 0) && (bedpeEntry_a.start1 == bedpeEntry_b.start1) && (bedpeEntry_a.end1 == bedpeEntry_b.end1) && (bedpeEntry_a.chrom2.compare(bedpeEntry_b.chrom2) <= 0))
        || ((bedpeEntry_a.chrom1.compare(bedpeEntry_b.chrom1) == 0) && (bedpeEntry_a.start1 == bedpeEntry_b.start1) && (bedpeEntry_a.end1 == bedpeEntry_b.end1) && (bedpeEntry_a.chrom2.compare(bedpeEntry_b.chrom2) == 0) && (bedpeEntry_a.start2 < bedpeEntry_b.start2 || (bedpeEntry_a.start2 == bedpeEntry_b.start2 && bedpeEntry_a.end2 < bedpeEntry_b.end2)));
}

inline bool operator>(BEDPE bedpeEntry_a, BEDPE bedpeEntry_b) {
    return bedpeEntry_b < bedpeEntry_a;
}

inline bool operator==(BEDPE bedpeEntry_a, BEDPE bedpeEntry_b) {
    return bedpeEntry_a.chrom1.compare(bedpeEntry_b.chrom1) == 0 && bedpeEntry_a.start1 == bedpeEntry_b.start1 && bedpeEntry_a.end1 == bedpeEntry_b.end1 &&
           bedpeEntry_a.chrom2.compare(bedpeEntry_b.chrom2) == 0 && bedpeEntry_a.start2 == bedpeEntry_b.start2 && bedpeEntry_a.end2 == bedpeEntry_b.end2;
}