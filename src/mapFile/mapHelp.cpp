/*
 * mapHelp.cpp
 *
 *  Created on: Apr 22, 2015
 *      Author: nek3d
 */

#include "CommonHelp.h"
#include "KeyListOps.h"

void map_help(void) {

    cerr << "\nTool:    bedtools map (aka mapBed)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Apply a function to a column from B intervals that overlap A." << endl << endl;

    cerr << "Usage:   " << "bedtools map" << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;

    cerr << "\t-c\t"             << "Specify columns from the B file to map onto intervals in A." << endl;
    cerr                         << "\t\tDefault: 5." << endl;
    cerr						<<  "\t\tMultiple columns can be specified in a comma-delimited list." << endl << endl;

    KeyListOpsHelp();


    cerr << "\t-f\t"             << "Minimum overlap required as a fraction of A." << endl;
    cerr                         << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                         << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-r\t"             << "Require that the fraction overlap be reciprocal for A and B." << endl;
    cerr                         << "\t\t- In other words, if -f is 0.90 and -r is used, this requires" << endl;
    cerr                         << "\t\t  that B overlap 90% of A and A _also_ overlaps 90% of B." << endl << endl;

    cerr << "\t-s\t"             << "Require same strandedness.  That is, only report hits in B" << endl;
    cerr                         << "\t\tthat overlap A on the _same_ strand." << endl;
    cerr                         << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-S\t"             << "Require different strandedness.  That is, only report hits in B" << endl;
    cerr                         << "\t\tthat overlap A on the _opposite_ strand." << endl;
    cerr                         << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl << endl;

    cerr << "\t-g\t"             << "Provide a genome file to enforce consistent chromosome sort order" << endl;
    cerr                         <<"\t\tacross input files." << endl << endl;

    cerr << "\t-null\t"          << "The value to print if no overlaps are found for an A interval." << endl;
    cerr                         << "\t\t- Default - \".\"" << endl << endl;

    cerr << "\t-header\t"        << "Print the header from the A file prior to results." << endl << endl;

    cerr << "\t-prec\t"   		 << "Sets the decimal precision for output (Default: 5)" << endl << endl;

    cerr << "\t-nonamecheck\t"       << "For sorted data, don't throw an error if the file has different naming conventions" << endl;
    cerr							<< "\t\t\tfor the same chromosome. ex. \"chr1\" vs \"chr01\"." << endl << endl;

    CommonHelp();

    cerr << "Notes: " << endl;
    cerr << "\t(1) Both input files must be sorted by chrom, then start." << endl << endl;

    // end the program here
    exit(1);

}



