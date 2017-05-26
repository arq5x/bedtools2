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
    cerr << "Summary: Returns the depth and breadth of coverage of features from B" << endl;
    cerr << "\t on the intervals in A." << endl << endl;

    cerr << "Usage:   " << "bedtools coverage" << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-hist\t"         << "Report a histogram of coverage for each feature in A" << endl;
    cerr                        << "\t\tas well as a summary histogram for _all_ features in A." << endl << endl;
    cerr                        << "\t\tOutput (tab delimited) after each feature in A:" << endl;
    cerr                        << "\t\t  1) depth\n\t\t  2) # bases at depth\n\t\t  3) size of A\n\t\t  4) % of A at depth" << endl << endl;

    cerr << "\t-d\t"            << "Report the depth at each position in each A feature." << endl;
    cerr                        << "\t\tPositions reported are one based.  Each position" << endl;
    cerr                        << "\t\tand depth follow the complete A feature." << endl << endl;
    
    cerr << "\t-counts\t"       << "Only report the count of overlaps, don't compute fraction, etc." << endl << endl;

    cerr << "\t-mean\t"       	<< "Report the mean depth of all positions in each A feature." << endl << endl;


    IntersectCommonHelp();
    sortedHelp();
    allToolsCommonHelp();

    cerr << "Default Output:  " << endl;
    cerr << "\t" << " After each entry in A, reports: " << endl;
    cerr << "\t   1) The number of features in B that overlapped the A interval." << endl;
    cerr << "\t   2) The number of bases in A that had non-zero coverage." << endl;
    cerr << "\t   3) The length of the entry in A." << endl;
    cerr << "\t   4) The fraction of bases in A that had non-zero coverage." << endl << endl;
}
