/*****************************************************************************
  closestMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "closestBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools closest"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void closest_help(void);

int closest_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedAFile;
    string bedBFile;
    string tieMode = "all";
    string strandedDistMode = "";

    bool haveBedA       = false;
    bool haveBedB       = false;
    bool haveTieMode    = false;
    bool sameStrand     = false;
    bool diffStrand     = false;
    bool ignoreOverlaps = false;
    bool ignoreUpstream = false;
    bool ignoreDownstream = false;
    bool reportDistance = false;
    bool signDistance   = false;
    bool haveStrandedDistMode = false;
    bool printHeader        = false;


    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if( (PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) closest_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-a", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveBedA = true;
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
        else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
            sameStrand = true;
        }
        else if (PARAMETER_CHECK("-S", 2, parameterLength)) {
            diffStrand = true;
        }
        else if (PARAMETER_CHECK("-d", 2, parameterLength)) {
            reportDistance = true;
        }
        else if (PARAMETER_CHECK("-D", 2, parameterLength)) {
            if ((i+1) < argc) {
                reportDistance = true;
                signDistance = true;
                haveStrandedDistMode = true;
                strandedDistMode = argv[i + 1];
                i++;
            }
        }
        else if (PARAMETER_CHECK("-io", 3, parameterLength)) {
            ignoreOverlaps = true;
        }
        else if (PARAMETER_CHECK("-iu", 3, parameterLength)) {
            ignoreUpstream = true;
        }
        else if (PARAMETER_CHECK("-id", 3, parameterLength)) {
            ignoreDownstream = true;
        }
        else if (PARAMETER_CHECK("-t", 2, parameterLength)) {
            if ((i+1) < argc) {
                haveTieMode = true;
                tieMode = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-header", 7, parameterLength)) {
            printHeader = true;
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

    if (haveTieMode && (tieMode != "all") && (tieMode != "first")
                    && (tieMode != "last")) {
        cerr << endl << "*****" << endl << "*****ERROR: Request \"all\" or \"first\" or \"last\" for Tie Mode (-t)" << endl << "*****" << endl;
        showHelp = true;
    }
    
    if (haveStrandedDistMode && (strandedDistMode != "a") && (strandedDistMode != "b")
                             && (strandedDistMode != "ref")) {
        cerr << endl << "*****" << endl << "*****ERROR: Request \"a\" or \"b\" or \"ref\" for Stranded Distance Mode (-D)" << endl << "*****" << endl;
        showHelp = true;
    }

    if (sameStrand && diffStrand) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -s OR -S, not both." << endl << "*****" << endl;
        showHelp = true;
    }
    
    if (ignoreUpstream && ignoreDownstream) {
        cerr << endl << "*****" << endl << "*****ERROR: Request either -iu OR -id, not both." << endl << "*****" << endl;
        showHelp = true;
    }
    
    if ((ignoreUpstream || ignoreDownstream) && ! haveStrandedDistMode) {
        cerr << endl << "*****" << endl << "*****ERROR: When requesting -iu or -id, you also need to specify -D." << endl << "*****" << endl;
        showHelp = true;
    }
    
    
    if (!showHelp) {
        BedClosest *bc = new BedClosest(bedAFile, bedBFile, sameStrand, 
                                        diffStrand, tieMode, reportDistance, 
                                        signDistance, strandedDistMode, ignoreOverlaps,
                                        ignoreUpstream, ignoreDownstream, printHeader);
        delete bc;
    }
    else {
        closest_help();
    }
    return 0;
}

void closest_help(void) {
    
    cerr << "\nTool:    bedtools closest (aka closestBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: For each feature in A, finds the closest " << endl;
    cerr << "\t feature (upstream or downstream) in B." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-s\t"            << "Req. same strandedness.  That is, find the closest feature in" << endl;
    cerr                        << "\t\tB that overlaps A on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-S\t"            << "Req. opposite strandedness.  That is, find the closest feature" << endl;
    cerr                        << "\t\tin B that overlaps A on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-d\t"            << "In addition to the closest feature in B, " << endl;
    cerr                        << "\t\treport its distance to A as an extra column." << endl;
    cerr                        << "\t\t- The reported distance for overlapping features will be 0." << endl << endl;
    
    cerr << "\t-D\t"            << "Like -d, report the closest feature in B, and its distance to A" << endl;
    cerr                        << "\t\tas an extra column. Unlike -d, use negative distances to report" << endl;
    cerr                        << "\t\tupstream features." << endl;
    cerr                        << "\t\tThe options for defining which orientation is \"upstream\" are:" << endl;
    cerr                        << "\t\t- \"ref\"   Report distance with respect to the reference genome. " << endl;
    cerr                        << "\t\t            B features with a lower (start, stop) are upstream" << endl;
    cerr                        << "\t\t- \"a\"     Report distance with respect to A." << endl;
    cerr                        << "\t\t            When A is on the - strand, \"upstream\" means B has a" << endl;
    cerr                        << "\t\t            higher (start,stop)." << endl;    
    cerr                        << "\t\t- \"b\"     Report distance with respect to B." << endl;
    cerr                        << "\t\t            When B is on the - strand, \"upstream\" means A has a" << endl;
    cerr                        << "\t\t            higher (start,stop)." << endl << endl;

    cerr << "\t-io\t"           << "Ignore features in B that overlap A.  That is, we want close," << endl;
    cerr                        << "\t\tyet not touching features only." << endl << endl;
    
    cerr << "\t-iu\t"           << "Ignore features in B that are upstream of features in A." << endl;
    cerr                        << "\t\tThis option requires -D and follows its orientation" << endl;
    cerr                        << "\t\trules for determining what is \"upstream\"." << endl;
    cerr << "\t-id\t"           << "Ignore features in B that are upstream of features in A." << endl;
    cerr                        << "\t\tThis option requires -D and follows its orientation" << endl;
    cerr                        << "\t\trules for determining what is \"downstream\"." << endl;

    cerr << "\t-t\t"            << "How ties for closest feature are handled.  This occurs when two" << endl;
    cerr                        << "\t\tfeatures in B have exactly the same \"closeness\" with A." << endl;
    cerr                        << "\t\tBy default, all such features in B are reported." << endl;
    cerr                        << "\t\tHere are all the options:" << endl;
    cerr                        << "\t\t- \"all\"    Report all ties (default)." << endl;
    cerr                        << "\t\t- \"first\"  Report the first tie that occurred in the B file." << endl;
    cerr                        << "\t\t- \"last\"   Report the last tie that occurred in the B file." << endl << endl;
    
    cerr << "\t-header\t"       << "Print the header from the A file prior to results." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\tReports \"none\" for chrom and \"-1\" for all other fields when a feature" << endl;
    cerr << "\tis not found in B on the same chromosome as the feature in A." << endl;
    cerr << "\tE.g. none\t-1\t-1" << endl << endl;

    // end the program here
    exit(1);
}
