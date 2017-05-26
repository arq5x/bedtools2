/*
 * sampleMain.cpp
 *
 *  Created on: Nov 18, 2013
 *      Author: nek3d
 */

#include "CommonHelp.h"

void sample_help(void) {

    cerr << "\nTool:    bedtools sample (aka sampleFile)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Take sample of input file(s) using reservoir sampling algorithm." << endl << endl;

    cerr << "Usage:   " << "bedtools sample" << " [OPTIONS] -i <bed/gff/vcf/bam>" << endl << endl;

    cerr << "WARNING:\tThe current sample algorithm will hold all requested sample records in memory prior to output." << endl;
    cerr << "\t\tThe user must ensure that there is adequate memory for this." << endl << endl;
    cerr << "Options: " << endl;

    cerr << "\t-n\t"                << "The number of records to generate." << endl;
    cerr                            << "\t\t- Default = 1,000,000." << endl;
    cerr                            << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-seed\t"             << "Supply an integer seed for the shuffling." << endl;
    cerr                            << "\t\t- By default, the seed is chosen automatically." << endl;
    cerr                            << "\t\t- (INTEGER)" << endl << endl;


    cerr << "\t-ubam\t"         << "Write uncompressed BAM output. Default writes compressed BAM." << endl << endl;

    cerr << "\t-s\t"            << "Require same strandedness.  That is, only give records" << endl;
    cerr                        << "\t\tthat have the same strand. Use '-s forward' or '-s reverse'" << endl;
    cerr						<< "\t\tfor forward or reverse strand records, respectively." << endl;
    cerr                        << "\t\t- By default, records are reported without respect to strand." << endl << endl;

    cerr << "\t-header\t"       << "Print the header from the input file prior to results." << endl << endl;

    allToolsCommonHelp();

    cerr << "Notes: " << endl;
}


