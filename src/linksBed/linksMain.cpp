/*****************************************************************************
  linksBedMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "linksBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools links"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void links_help(void);

int links_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;

    // input files
    string bedFile = "stdin";
    bool haveBed   = true;

    /* Defaults for everyone else */
    string org = "human";
    string db = "hg18";
    string base = "http://genome.ucsc.edu";

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) links_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedFile = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-base", 5, parameterLength)) {
            if ((i+1) < argc) {
                base = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-org", 4, parameterLength)) {
            if ((i+1) < argc) {
                org = argv[i + 1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-db", 3, parameterLength)) {
            if ((i+1) < argc) {
                db = argv[i + 1];
                i++;
            }
        }
        else {
            cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
            showHelp = true;
        }
    }

    // make sure we have both input files
    if (!haveBed) {
        cerr << endl << "*****" << endl << "*****ERROR: Need -i BED file. " << endl << "*****" << endl;
        showHelp = true;
    }

    if (!showHelp) {
        BedLinks *bl = new BedLinks(bedFile, base, org, db);
        delete bl;
    }
    else {
        links_help();
    }
    return 0;
}

void links_help(void) {

    cerr << "\nTool:    bedtools links (aka linksBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Creates HTML links to an UCSC Genome Browser from a feature file." << endl << endl;
    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf> > out.html" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-base\t" << "The browser basename.  Default: http://genome.ucsc.edu " << endl;
    cerr << "\t-org\t"  << "The organism. Default: human" << endl;
    cerr << "\t-db\t"   << "The build.  Default: hg18" << endl << endl;

    cerr << "Example: " << endl;
    cerr << "\t" << "By default, the links created will point to human (hg18) UCSC browser." << endl;
    cerr <<         "\tIf you have a local mirror, you can override this behavior by supplying" << endl;
    cerr <<         "\tthe -base, -org, and -db options."  << endl << endl;
    cerr << "\t" << "For example, if the URL of your local mirror for mouse MM9 is called: " << endl;
    cerr <<         "\thttp://mymirror.myuniversity.edu, then you would use the following:" << endl;
    cerr <<         "\t" << "-base http://mymirror.myuniversity.edu" << endl;
    cerr <<         "\t" << "-org mouse" << endl;
    cerr <<         "\t" << "-db mm9" << endl;


    exit(1);
}
