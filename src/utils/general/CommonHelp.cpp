/*
 * CommonHelp.cpp
 *
 *  Created on: Jul 2, 2014
 *      Author: nek3d
 */
#include "CommonHelp.h"

void IntersectOutputHelp() {
	cerr << "\t-wa\t"           << "Write the original entry in A for each overlap." << endl << endl;

	cerr << "\t-wb\t"           << "Write the original entry in B for each overlap." << endl;
	cerr                        << "\t\t- Useful for knowing _what_ A overlaps. Restricted by -f and -r." << endl << endl;

	cerr << "\t-loj\t"          << "Perform a \"left outer join\". That is, for each feature in A" << endl;
	cerr                        << "\t\treport each overlap with B.  If no overlaps are found, " << endl;
	cerr                        << "\t\treport a NULL feature for B." << endl << endl;

	cerr << "\t-wo\t"           << "Write the original A and B entries plus the number of base" << endl;
	cerr                        << "\t\tpairs of overlap between the two features." << endl;
	cerr                        << "\t\t- Overlaps restricted by -f and -r." << endl;
	cerr                        << "\t\t  Only A features with overlap are reported." << endl << endl;

	cerr << "\t-wao\t"          << "Write the original A and B entries plus the number of base" << endl;
	cerr                        << "\t\tpairs of overlap between the two features." << endl;
	cerr                        << "\t\t- Overlapping features restricted by -f and -r." << endl;
	cerr                        << "\t\t  However, A features w/o overlap are also reported" << endl;
	cerr                        << "\t\t  with a NULL B feature and overlap = 0." << endl << endl;

	cerr << "\t-u\t"            << "Write the original A entry _once_ if _any_ overlaps found in B." << endl;
	cerr                        << "\t\t- In other words, just report the fact >=1 hit was found." << endl;
	cerr                        << "\t\t- Overlaps restricted by -f and -r." << endl << endl;

    cerr << "\t-c\t"            << "For each entry in A, report the number of overlaps with B." << endl;
    cerr                        << "\t\t- Reports 0 for A entries that have no overlap with B." << endl;
    cerr                        << "\t\t- Overlaps restricted by -f, -F, -r, and -s." << endl << endl;

    cerr << "\t-C\t"            << "For each entry in A, separately report the number of" << endl;
    cerr                        << "\t\t- overlaps with each B file on a distinct line." << endl;
    cerr                        << "\t\t- Reports 0 for A entries that have no overlap with B." << endl;
    cerr                        << "\t\t- Overlaps restricted by -f, -F, -r, and -s." << endl << endl;

    cerr << "\t-v\t"            << "Only report those entries in A that have _no overlaps_ with B." << endl;
    cerr                        << "\t\t- Similar to \"grep -v\" (an homage)." << endl << endl;

    cerr << "\t-ubam\t"         << "Write uncompressed BAM output. Default writes compressed BAM." << endl << endl;

}

void sortedHelp() {
    cerr << "\t-sorted\t"       << "Use the \"chromsweep\" algorithm for sorted (-k1,1 -k2,2n) input." << endl << endl;
}


void IntersectCommonHelp() {
    cerr << "\t-s\t"            << "Require same strandedness.  That is, only report hits in B" << endl;
    cerr                        << "\t\tthat overlap A on the _same_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-S\t"            << "Require different strandedness.  That is, only report hits in B" << endl;
    cerr                        << "\t\tthat overlap A on the _opposite_ strand." << endl;
    cerr                        << "\t\t- By default, overlaps are reported without respect to strand." << endl << endl;

    cerr << "\t-f\t"            << "Minimum overlap required as a fraction of A." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

	cerr << "\t-F\t"            << "Minimum overlap required as a fraction of B." << endl;
    cerr                        << "\t\t- Default is 1E-9 (i.e., 1bp)." << endl;
    cerr                        << "\t\t- FLOAT (e.g. 0.50)" << endl << endl;

    cerr << "\t-r\t"            << "Require that the fraction overlap be reciprocal for A AND B." << endl;
    cerr                        << "\t\t- In other words, if -f is 0.90 and -r is used, this requires" << endl;
    cerr                        << "\t\t  that B overlap 90% of A and A _also_ overlaps 90% of B." << endl << endl;

	cerr << "\t-e\t"            << "Require that the minimum fraction be satisfied for A OR B." << endl;
    cerr                        << "\t\t- In other words, if -e is used with -f 0.90 and -F 0.10 this requires" << endl;
    cerr                        << "\t\t  that either 90% of A is covered OR 10% of  B is covered." << endl;
    cerr                        << "\t\t  Without -e, both fractions would have to be satisfied." << endl << endl;

    cerr << "\t-split\t"        << "Treat \"split\" BAM or BED12 entries as distinct BED intervals." << endl << endl;

    cerr << "\t-g\t"       		<< "Provide a genome file to enforce consistent chromosome sort order" << endl;
    cerr 						<<"\t\tacross input files. Only applies when used with -sorted option." << endl << endl;

    cerr << "\t-nonamecheck\t"       << "For sorted data, don't throw an error if the file has different naming conventions" << endl;
    cerr							<< "\t\t\tfor the same chromosome. ex. \"chr1\" vs \"chr01\"." << endl << endl;

}

void multiDbOutputHelp() {
    cerr << "\t-names\t"       << "When using multiple databases, provide an alias for each that" << endl;
    cerr						<<"\t\twill appear instead of a fileId when also printing the DB record." << endl << endl;

    cerr << "\t-filenames"       << "\tWhen using multiple databases, show each complete filename" << endl;
    cerr						<<"\t\t\tinstead of a fileId when also printing the DB record." << endl << endl;

    cerr << "\t-sortout\t"       << "When using multiple databases, sort the output DB hits" << endl;
    cerr						<< "\t\t\tfor each record." << endl << endl;
}


void allToolsCommonHelp() {
	cerr << "\t-bed\t"          << "If using BAM input, write output as BED." << endl << endl;

	cerr << "\t-header\t"       << "Print the header from the A file prior to results." << endl << endl;

	cerr << "\t-nobuf\t"       << "Disable buffered output. Using this option will cause each line"<< endl;
	cerr 						<<"\t\tof output to be printed as it is generated, rather than saved" << endl;
	cerr 						<<"\t\tin a buffer. This will make printing large output files " << endl;

	cerr 						<<"\t\tnoticeably slower, but can be useful in conjunction with" << endl;
	cerr 						<<"\t\tother software tools and scripts that need to process one" << endl;
	cerr 						<<"\t\tline of bedtools output at a time." << endl << endl;

	cerr << "\t-iobuf\t"            << "Specify amount of memory to use for input buffer." << endl;
	cerr << "\t\t" <<					"Takes an integer argument. Optional suffixes K/M/G supported." << endl;
	cerr << "\t\t" 					<< "Note: currently has no effect with compressed files." << endl << endl;
}
