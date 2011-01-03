/*****************************************************************************
  intersectMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "intersectBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "intersectBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;

    // input arguments
    float overlapFraction = 1E-9;

    bool haveBedA           = false;
    bool haveBedB           = false;
    bool noHit              = false;
    bool anyHit             = false;
    bool writeA             = false;
    bool writeB             = false;
    bool writeCount         = false;
    bool writeOverlap       = false;
    bool writeAllOverlap    = false;
    bool haveFraction       = false;
    bool reciprocalFraction = false;
    bool forceStrand        = false;
    bool obeySplits         = false;
    bool inputIsBam         = false;
    bool outputIsBam        = true;
    bool uncompressedBam    = false;
    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

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

        if(PARAMETER_CHECK("-a", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBedA = true;
                outputIsBam = false;
                bedAFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-abam", 5, parameterLength)) {
            if ((i+1) < argc) {
                haveBedA = true;
                inputIsBam = true;
                bedAFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBedB = true;
                bedBFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-bed", 4, parameterLength)) {
            outputIsBam = false;
        }
        else if(PARAMETER_CHECK("-u", 2, parameterLength)) {
            anyHit = true;
        }
        else if(PARAMETER_CHECK("-f", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveFraction = true;
                overlapFraction = atof(argv[i + 1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-wa", 3, parameterLength)) {
            writeA = true;
        }
        else if(PARAMETER_CHECK("-wb", 3, parameterLength)) {
            writeB = true;
        }
        else if(PARAMETER_CHECK("-wo", 3, parameterLength)) {
            writeOverlap = true;
        }
        else if(PARAMETER_CHECK("-wao", 4, parameterLength)) {
            writeAllOverlap = true;
            writeOverlap = true;
        }
        else if(PARAMETER_CHECK("-c", 2, parameterLength)) {
            writeCount = true;
        }
        else if(PARAMETER_CHECK("-r", 2, parameterLength)) {
            reciprocalFraction = true;
        }
        else if (PARAMETER_CHECK("-v", 2, parameterLength)) {
            noHit = true;
        }
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            forceStrand = true;
        }
        else if (PARAMETER_CHECK("-split", 6, parameterLength)) {
            obeySplits = true;
        }
        else if(PARAMETER_CHECK("-ubam", 5, parameterLength)) {
            uncompressedBam = true;
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBedA || !haveBedB) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -a and -b files. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (anyHit && noHit) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -v, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (writeB && writeCount) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -wb OR -c, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (writeCount && writeOverlap) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -wb OR -wo, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (writeA && writeOverlap) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -wa OR -wo, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (writeB && writeOverlap) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -wb OR -wo, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (reciprocalFraction && !haveFraction) {
        cerr << endl << "*****" << endl << "*****ERROR: If using -r, you need to define -f." << endl << "*****" << endl;
        showHelp = true;
    }

    if (anyHit && writeCount) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -c, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (anyHit && writeB) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -wb, not both." << endl << "*****" << endl;
        showHelp = true;
    }

    if (anyHit && writeOverlap) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -wo, not both." << endl << "*****" << endl;
        showHelp = true;
    }


    if (!showHelp) {

        BedIntersect *bi = new BedIntersect(bedAFile, bedBFile, anyHit, writeA, writeB, writeOverlap,
                                            writeAllOverlap, overlapFraction, noHit, writeCount, forceStrand,
                                            reciprocalFraction, obeySplits, inputIsBam, outputIsBam, uncompressedBam);
        delete bi;
        return 0;
    }
    else {
        ShowHelp();
    }
}

void ShowHelp(void) {

    cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;

    cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

    cerr << "Summary: Report overlaps between two feature files." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-abam\t"         << "The A input file is in BAM format.  Output will be BAM as well." << endl << endl;

    cerr << "\t-ubam\t"         << "Write uncompressed BAM output. Default is to write compressed BAM." << endl << endl;

    cerr << "\t-bed\t"          << "When using BAM input (-abam), write output as BED. The default" << endl;
    cerr                        << "\t\tis to write output in BAM when using -abam." << endl << endl;

    cerr << "\t-wa\t"           << "Write the original entry in A for each overlap." << endl << endl;

    cerr << "\t-wb\t"           << "Write the original entry in B for each overlap." << endl;
    cerr                        << "\t\t- Useful for knowing _what_ A overlaps. Restricted by -f and -r." << endl << endl;

    cerr << "\t-wo\t"           << "Write the original A and B entries plus the number of base" << endl;
    cerr                        << "\t\tpairs of overlap between the two features." << endl;
    cerr                        << "\t\t- Overlaps restricted by -f and -r." << endl;
    cerr                        << "\t\t  Only A features with overlap are reported." << endl << endl;

    cerr << "\t-wao\t"          << "Write the original A and B entries plus the number of base" << endl;
    cerr                        << "\t\tpairs of overlap between the two features." << endl;
    cerr                        << "\t\t- Overlapping features restricted by -f and -r." << endl;
    cerr                        << "\t\t  However, A features w/o overlap are also reported" << endl;
    cerr                        << "\t\t  with a NULL B feature and overlap = 0." << endl << endl;

    cerr << "\t-u\t"            << "Write the original A entry _once_ if _any_ overlaps found in B." << endl;
    cerr                        << "\t\t- In other words, just report the fact >=1 hit was found." << endl;
    cerr                        << "\t\t- Overlaps restricted by -f and -r." << endl << endl;

    cerr << "\t-c\t"            << "For each entry in A, report the number of overlaps with B." << endl;
    cerr                        << "\t\t- Reports 0 for A entries that have no overlap with B." << endl;
    cerr                        << "\t\t- Overlaps restricted by -f and -r." << endl << endl;

    cerr << "\t-v\t"            << "Only report those entries in A that have _no overlaps_ with B." << endl;
    cerr                        << "\t\t- Similar to \"grep -v\" (an homage)." << endl << endl;

    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of A." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-r\t"            << "Require that the fraction overlap be reciprocal for A and B." << endl;
    cerr                        << "\t\t- In other words, if -f is 0.90 and -r is used, this requires" << endl;
    cerr                        << "\t\t  that B overlap 90% of A and A _also_ overlaps 90% of B." << endl << endl;

    cerr << "\t-s\t"            << "Force strandedness.  That is, only report hits in B that" << endl;
    cerr                        << "\t\toverlap A on the same strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl << endl;


    // end the program here
    exit(1);

}
