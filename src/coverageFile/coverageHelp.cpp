/*****************************************************************************
  coverageMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "CommonHelp.h"


void coverage_help(void) {

    cerr << "\nTool:    bedtools coverage (aka coverageBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Returns the depth and breadth of coverage of features from A" << endl;
    cerr << "\t on the intervals in B." << endl << endl;

    cerr << "Usage:   " << "bedtools coverage" << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-abam\t"         << "The A input file is in BAM format. Replaces -a." << endl << endl;

    cerr << "\t-s\t"            << "Require same strandedness.  That is, only counts hits in A that" << endl;
    cerr                        << "\t\toverlap B on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are counted without respect to strand." << endl << endl;

    cerr << "\t-S\t"            << "Require different strandedness.  That is, only report hits in A" << endl;
    cerr                        << "\t\tthat overlap B on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are counted without respect to strand." << endl << endl;

    cerr << "\t-hist\t"         << "Report a histogram of coverage for each feature in B" << endl;
    cerr                        << "\t\tas well as a summary histogram for _all_ features in B." << endl << endl;
    cerr                        << "\t\tOutput (tab delimited) after each feature in B:" << endl;
    cerr                        << "\t\t  1) depth\n\t\t  2) # bases at depth\n\t\t  3) size of B\n\t\t  4) % of B at depth" << endl << endl;

    cerr << "\t-d\t"            << "Report the depth at each position in each B feature." << endl;
    cerr                        << "\t\tPositions reported are one based.  Each position" << endl;
    cerr                        << "\t\tand depth follow the complete B feature." << endl << endl;
    
    cerr << "\t-counts\t"       << "Only report the count of overlaps, don't compute fraction, etc." << endl << endl;

    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl;
    cerr                        << "\t\twhen computing coverage." << endl;
    cerr                        << "\t\tFor BAM files, this uses the CIGAR \"N\" and \"D\" operations " << endl;
    cerr                        << "\t\tto infer the blocks for computing coverage." << endl;
    cerr                        << "\t\tFor BED12 files, this uses the BlockCount, BlockStarts," << endl;
    cerr                        << "\t\tand BlockEnds fields (i.e., columns 10,11,12)." << endl << endl;


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

    CommonHelp();

    cerr << "Default Output:  " << endl;
    cerr << "\t" << " After each entry in B, reports: " << endl;
    cerr << "\t   1) The number of features in A that overlapped the B interval." << endl;
    cerr << "\t   2) The number of bases in B that had non-zero coverage." << endl;
    cerr << "\t   3) The length of the entry in B." << endl;
    cerr << "\t   4) The fraction of bases in B that had non-zero coverage." << endl << endl;

    exit(1);
}
