/*****************************************************************************
  jaccardMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "CommonHelp.h"

void jaccard_help(void) {

    cerr << "\nTool:    bedtools jaccard (aka jaccard)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Calculate Jaccard statistic b/w two feature files."
         << endl
         << "\t Jaccard is the length of the intersection over the union."
         << endl
         << "\t Values range from 0 (no intersection) to 1 (self intersection)."
         << endl << endl;

    cerr << "Usage:   " << "bedtools jaccard" << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;


    IntersectCommonHelp();
    allToolsCommonHelp();

    cerr << "Notes: " << endl;
    cerr << "\t(1) Input files must be sorted by chrom, then start position."
         << endl << endl;
}
