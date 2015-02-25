/*
 * subtractMain.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */


#include "subtractFile.h"
#include "ContextSubtract.h"
#include "CommonHelp.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools subtract"

void subtract_help(void);

int subtract_main(int argc, char* argv[]) {

    ContextSubtract *context = new ContextSubtract();
    if (!context->parseCmdArgs(argc, argv, 1) || context->getShowHelp() || !context->isValidState()) {
    	if (!context->getErrorMsg().empty()) {
    		cerr << context->getErrorMsg() << endl;
    	}
    	subtract_help();
    	delete context;
    	return 1;
    }
	SubtractFile *subtractFile = new SubtractFile(context);

	bool retVal = subtractFile->subtractFiles();
	delete subtractFile;
	delete context;
	return retVal ? 0 : 1;
}

void subtract_help(void) {
    cerr << "\nTool:    bedtools subtract (aka subtractBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Removes the portion(s) of an interval that is overlapped" << endl;
    cerr << "\t by another feature(s)." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of A." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- (FLOAT) (e.g. 0.50)" << endl << endl;

    cerr << "\t-s\t"            << "Require same strandedness.  That is, only subtract hits in B" << endl;
    cerr                        << "\t\tthat overlap A on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are subtracted without respect to strand." << endl << endl;

    cerr << "\t-S\t"            << "Force strandedness.  That is, only subtract hits in B that" << endl;
    cerr                        << "\t\toverlap A on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are subtracted without respect to strand." << endl << endl;

    cerr << "\t-A\t"            << "Remove entire feature if any overlap.  That is, by default," << endl;
    cerr                        << "\t\tonly subtract the portion of A that overlaps B. Here, if" << endl;
    cerr                        << "\t\tany overlap is found (or -f amount), the entire feature is removed." << endl << endl;

    cerr << "\t-N\t"            << "Same as -A except when used with -f, the amount is the sum" << endl;
    cerr                        << "\t\tof all features (not any single feature)." << endl << endl;


    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl << endl;

    cerr << "\t-sorted\t"       << "Use the \"chromsweep\" algorithm for sorted (-k1,1 -k2,2n) input." << endl << endl;

    cerr << "\t-g\t"       		<< "Provide a genome file to enforce consistent chromosome sort order" << endl;
    cerr 						<<"\t\tacross input files. Only applies when used with -sorted option." << endl << endl;

    cerr << "\t-header\t"       << "Print the header from the A file prior to results." << endl << endl;

    cerr << "\t-nobuf\t"       << "Disable buffered output. Using this option will cause each line"<< endl;
    cerr 						<<"\t\tof output to be printed as it is generated, rather than saved" << endl;
    cerr 						<<"\t\tin a buffer. This will make printing large output files " << endl;

    cerr 						<<"\t\tnoticeably slower, but can be useful in conjunction with" << endl;
    cerr 						<<"\t\tother software tools and scripts that need to process one" << endl;
    cerr 						<<"\t\tline of bedtools output at a time." << endl << endl;

    cerr << "\t-names\t"       << "When using multiple databases, provide an alias for each that" << endl;
    cerr						<<"\t\twill appear instead of a fileId when also printing the DB record." << endl << endl;

    cerr << "\t-filenames"       << "\tWhen using multiple databases, show each complete filename" << endl;
    cerr						<<"\t\t\tinstead of a fileId when also printing the DB record." << endl << endl;

    cerr << "\t-sortout\t"       << "When using multiple databases, sort the output DB hits" << endl;
    cerr						<< "\t\t\tfor each record." << endl << endl;

    cerr << "\t-nonamecheck\t"       << "For sorted data, don't throw an error if the file has different naming conventions" << endl;
    cerr							<< "\t\t\tfor the same chromosome. ex. \"chr1\" vs \"chr01\"." << endl << endl;

    // end the program here
    exit(1);

}
