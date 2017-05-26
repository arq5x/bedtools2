#include "CommonHelp.h"

void fisher_help(void) {

    cerr << "\nTool:    bedtools fisher (aka fisher)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Calculate Fisher statistic b/w two feature files."
         << endl << endl;

    cerr << "Usage:   " << "bedtools fisher" << " [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf> -g <genome file>" << endl << endl;

    cerr << "Options: " << endl;


    cerr << "\t-m\t"            << "Merge overlapping intervals before" << endl;
    cerr                        << "\t\t- looking at overlap." << endl << endl;

    IntersectCommonHelp();
    allToolsCommonHelp();


    cerr << "Notes: " << endl;
    cerr << "\t(1) Input files must be sorted by chrom, then start position."
         << endl << endl;
}
