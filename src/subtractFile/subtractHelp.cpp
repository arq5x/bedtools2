/*
 * subtractMain.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: nek3d
 */


#include "CommonHelp.h"

void subtract_help(void) {
    cerr << "\nTool:    bedtools subtract (aka subtractBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Removes the portion(s) of an interval that is overlapped" << endl;
    cerr << "\t by another feature(s)." << endl << endl;

    cerr << "Usage:   " << "bedtools subtract" << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-A\t"            << "Remove entire feature if any overlap.  That is, by default," << endl;
    cerr                        << "\t\tonly subtract the portion of A that overlaps B. Here, if" << endl;
    cerr                        << "\t\tany overlap is found (or -f amount), the entire feature is removed." << endl << endl;

    cerr << "\t-N\t"            << "Same as -A except when used with -f, the amount is the sum" << endl;
    cerr                        << "\t\tof all features (not any single feature)." << endl << endl;

    cerr << "\t-wb\t"           << "Write the original entry in B for each overlap." << endl;
    cerr                        << "\t\t- Useful for knowing _what_ A overlaps. Restricted by -f and -r." << endl << endl;

    cerr << "\t-wo\t"           << "Write the original A and B entries plus the number of base" << endl;
    cerr                        << "\t\tpairs of overlap between the two features." << endl;
    cerr                        << "\t\t- Overlaps restricted by -f and -r." << endl;
    cerr                        << "\t\t  Only A features with overlap are reported." << endl << endl;

    IntersectCommonHelp();
    sortedHelp();
    allToolsCommonHelp();
}
