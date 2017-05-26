/*****************************************************************************
  intersectHelp.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/

#include "CommonHelp.h"

void intersect_help(void) {

    cerr << "\nTool:    bedtools intersect (aka intersectBed)" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Report overlaps between two feature files." << endl << endl;

    cerr << "Usage:   " << "bedtools intersect" << " [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>" << endl << endl;

    cerr << "\t"				<< "Note: -b may be followed with multiple databases and/or " << endl;
    cerr << "\t"					"wildcard (*) character(s). " << endl;

    cerr << "Options: " << endl;

    // -abam is obsolete.
    // cerr << "\t-abam\t"         << "The A input file is in BAM format.  Output will be BAM as well." << endl << endl;


    IntersectOutputHelp();
    IntersectCommonHelp();

    sortedHelp();
    multiDbOutputHelp();
    allToolsCommonHelp();

    cerr << "Notes: " << endl;
    cerr << "\t(1) When a BAM file is used for the A file, the alignment is retained if overlaps exist," << endl;
    cerr << "\tand excluded if an overlap cannot be found.  If multiple overlaps exist, they are not" << endl;
    cerr << "\treported, as we are only testing for one or more overlaps." << endl << endl;

}
