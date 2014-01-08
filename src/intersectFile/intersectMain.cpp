/*****************************************************************************
  intersectMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
using namespace std;

#include "intersectFile.h"
#include "ContextIntersect.h"

// define our program name
#define PROGRAM_NAME "bedtools intersect"

void intersect_help(void);

int intersect_main(int argc, char* argv[]) {

    ContextIntersect *context = new ContextIntersect();
    if (!context->parseCmdArgs(argc, argv, 1) || context->getShowHelp() || !context->isValidState()) {
    	if (!context->getErrorMsg().empty()) {
    		cerr << context->getErrorMsg() << endl;
    	}
    	intersect_help();
    	delete context;
    	return 0;
    }
	FileIntersect *fileIntersect = new FileIntersect(context);

	bool retVal = fileIntersect->intersectFiles();
	delete fileIntersect;
	delete context;
	return retVal ? 0 : 1;
}

void intersect_help(void) {

    cerr << "\nTool:    bedtools intersect (aka intersectBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Report overlaps between two feature files." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-abam\t"         << "The A input file is in BAM format.  Output will be BAM as well." << endl << endl;

    cerr << "\t-ubam\t"         << "Write uncompressed BAM output. Default writes compressed BAM." << endl << endl;

    cerr << "\t-bed\t"          << "When using BAM input (-abam), write output as BED. The default" << endl;
    cerr                        << "\t\tis to write output in BAM when using -abam." << endl << endl;

    cerr << "\t-wa\t"           << "Write the original entry in A for each overlap." << endl << endl;

    cerr << "\t-wb\t"           << "Write the original entry in B for each overlap." << endl;
    cerr                        << "\t\t- Useful for knowing _what_ A overlaps. Restricted by -f and -r." << endl << endl;
    
    cerr << "\t-loj\t"          << "Perform a \"left outer join\". That is, for each feature in A" << endl;
    cerr                        << "\t\treport each overlap with B.  If no overlaps are found, " << endl;
    cerr                        << "\t\treport a NULL feature for B." << endl << endl;

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

    cerr << "\t-s\t"            << "Require same strandedness.  That is, only report hits in B" << endl;
    cerr                        << "\t\tthat overlap A on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-S\t"            << "Require different strandedness.  That is, only report hits in B" << endl;
    cerr                        << "\t\tthat overlap A on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

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

    cerr << "Notes: " << endl;
    cerr << "\t(1) When a BAM file is used for the A file, the alignment is retained if overlaps exist," << endl;
    cerr << "\tand exlcuded if an overlap cannot be found.  If multiple overlaps exist, they are not" << endl;
    cerr << "\treported, as we are only testing for one or more overlaps." << endl << endl;

    // end the program here
    exit(1);

}
