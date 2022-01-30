/*
 * qcHelp.cpp
 *
 *  Created on: Jan 2, 2019
 *      Author: Aaron Quinlan
 */

#include "CommonHelp.h"

void qc_help(void) {

    cerr << "\nTool:    bedtools qc" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Validate and report summary statistics of a file " << endl << endl;

    cerr << "Usage:   " << "bedtools qc" << " [OPTIONS] -i <bed/gff/vcf/bam> -g <genome>" << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\t(1)  The genome file should tab delimited and structured as follows:" << endl;
    cerr << "\t     <chromName><TAB><chromSize>" << endl << endl;
    cerr << "\tFor example, Human (hg19):" << endl;
    cerr << "\tchr1\t249250621" << endl;
    cerr << "\tchr2\t243199373" << endl;
    cerr << "\t..." << endl;
    cerr << "\tchr18_gl000207_random\t4262" << endl << endl;
}
