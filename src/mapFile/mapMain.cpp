/*****************************************************************************
  mapMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
using namespace std;

#include "mapFile.h"
#include "ContextMap.h"

// define our program name
#define PROGRAM_NAME "bedtools map"

void map_help(void);

int map_main(int argc, char* argv[]) {

    ContextMap *context = new ContextMap();
    if (!context->parseCmdArgs(argc, argv, 1) || context->getShowHelp() || !context->isValidState()) {
        if (!context->getErrorMsg().empty()) {
            cerr << context->getErrorMsg() << endl;
        }
        map_help();
        delete context;
        return 0;
    }
    FileMap *fileMap = new FileMap(context);

    bool retVal = fileMap->mapFiles();
    delete fileMap;
    delete context;
    return retVal ? 0 : 1;
}

void map_help(void) {

    cerr << "\nTool:    bedtools map (aka mapBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Apply a function to a column from B intervals that overlap A." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-c\t"             << "Specify columns from the B file to map onto intervals in A." << endl;
    cerr                         << "\t\tDefault: 5." << endl;
    cerr						<<  "\t\tMultiple columns can be specified in a comma-delimited list." << endl << endl;

    KeyListOpsHelp();


    cerr << "\t-f\t"             << "Minimum overlap required as a fraction of A." << endl;
    cerr                         << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                         << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;
                                 
    cerr << "\t-r\t"             << "Require that the fraction overlap be reciprocal for A and B." << endl;
    cerr                         << "\t\t- In other words, if -f is 0.90 and -r is used, this requires" << endl;
    cerr                         << "\t\t  that B overlap 90% of A and A _also_ overlaps 90% of B." << endl << endl;
                                 
    cerr << "\t-s\t"             << "Require same strandedness.  That is, only report hits in B" << endl;
    cerr                         << "\t\tthat overlap A on the _same_ strand." << endl;
    cerr                         << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;
                                 
    cerr << "\t-S\t"             << "Require different strandedness.  That is, only report hits in B" << endl;
    cerr                         << "\t\tthat overlap A on the _opposite_ strand." << endl;
    cerr                         << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl << endl;

    cerr << "\t-g\t"             << "Provide a genome file to enforce consistent chromosome sort order" << endl;
    cerr                         <<"\t\tacross input files." << endl << endl;

    cerr << "\t-null\t"          << "The value to print if no overlaps are found for an A interval." << endl;
    cerr                         << "\t\t- Default - \".\"" << endl << endl;

    cerr << "\t-header\t"        << "Print the header from the A file prior to results." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1) Both input files must be sorted by chrom, then start." << endl << endl;
    
    // end the program here
    exit(1);

}
