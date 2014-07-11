#include "Fisher.h"
#include "version.h"

using namespace std;

#define PROGRAM_NAME "bedtools fisher"

void fisher_help(void);

int fisher_main(int argc, char* argv[]) {

    ContextFisher *context = new ContextFisher();
    if (!context->parseCmdArgs(argc, argv, 1) || context->getShowHelp() || !context->isValidState()) {
    	if (!context->getErrorMsg().empty()) {
    		cerr << context->getErrorMsg() << endl;
    	}
    	fisher_help();
    	delete context;
    	return 0;
    }
	Fisher *fisher = new Fisher(context);

	bool retVal = fisher->calculate();
	delete fisher;
	delete context;
	return retVal ? 0 : 1;
}


void fisher_help(void) {

    cerr << "\nTool:    bedtools fisher (aka fisher)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Calculate Fisher statistic b/w two feature files."
         << endl
         << "\t Fisher is the length of the intersection over the union."
         << endl
         << "\t Values range from 0 (no intersection) to 1 (self intersection)."
         << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf> -g <genome>" << endl << endl;

    cerr << "Options: " << endl;


    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of A." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-r\t"            << "Require that the fraction overlap be reciprocal for A and B." << endl;
    cerr                        << "\t\t- In other words, if -f is 0.90 and -r is used, this requires" << endl;
    cerr                        << "\t\t  that B overlap 90% of A and A _also_ overlaps 90% of B." << endl << endl;

    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl << endl;

    cerr << "\t-s\t"           << "Force strandedness.  That is, only merge features" << endl;
    cerr                       << "\t\tthat are on the same strand." << endl;
    cerr                       << "\t\t- By default, merging is done without respect to strand." << endl << endl;

    cerr << "\t-S\t"		   << "Force merge for one specific strand only." << endl;
    cerr << "\t\t"             << "Follow with + or - to force merge from only" << endl;
    cerr << "\t\t"			   << "the forward or reverse strand, respectively." << endl;
    cerr << "\t\t"			   << "- By default, merging is done without respect to strand." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1) Input files must be sorted by chrom, then start position."
         << endl << endl;

    // end the program here
    exit(1);

}
