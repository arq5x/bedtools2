/*
 * spacingMain.cpp
 *
 *  Created on: Jan 15, 2015
 *      Author: Aaron Quinlan
 */

#include "CommonHelp.h"

void spacing_help(void) {

    cerr << "\nTool:    bedtools spacing" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Report (last col.) the gap lengths between intervals in a file." << endl << endl;

    cerr << "Usage:   " << "bedtools spacing" << " [OPTIONS] -i <bed/gff/vcf/bam>" << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  Input must be sorted by chrom,start (sort -k1,1 -k2,2n for BED)." << endl;
    cerr << "\t(2)  The 1st element for each chrom will have NULL distance. (\".\")." << endl;
    cerr << "\t(3)  Distance for overlapping intervals is -1 and 0 for adjacent intervals." << endl << endl;

    cerr << "Example: " << endl;
    cerr << "\t$ cat test.bed " << endl;
    cerr << "\tchr1    0   10 " << endl;
    cerr << "\tchr1    10  20 " << endl;
    cerr << "\tchr1    19  30 " << endl;
    cerr << "\tchr1    35  45 " << endl;
    cerr << "\tchr1    100 200 " << endl << endl;

    cerr << "\t$ bedtools spacing -i test.bed " << endl;
    cerr << "\tchr1    0   10  . " << endl;
    cerr << "\tchr1    10  20  0 " << endl;
    cerr << "\tchr1    19  30  -1 " << endl;
    cerr << "\tchr1    35  45  5 " << endl;
    cerr << "\tchr1    100 200 55 " << endl << endl;

    allToolsCommonHelp();
}
