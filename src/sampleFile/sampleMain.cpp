/*
 * sampleMain.cpp
 *
 *  Created on: Nov 18, 2013
 *      Author: nek3d
 */

#include <iostream>
#include "ContextSample.h"
#include "SampleFile.h"

#define PROGRAM_NAME "bedtools sample"

void sample_help(void);

int sample_main(int argc, char **argv)
{
    ContextSample *context = new ContextSample();
    if (!context->parseCmdArgs(argc, argv, 1) || context->getShowHelp() || !context->isValidState()) {
    	if (!context->getErrorMsg().empty()) {
    		cerr << context->getErrorMsg() << endl;
    	}
    	sample_help();
    	delete context;
    	return 0;
    }
    SampleFile *sampleFile = new SampleFile(context);

	bool retVal = sampleFile->takeSample();
	delete sampleFile;
	delete context;
	return retVal ? 0 : 1;

}

void sample_help(void) {

    cerr << "\nTool:    bedtools sample (aka sampleFile)" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Take sample of input file(s) using reservoir sampling algorithm." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf/bam>" << endl << endl;

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

    cerr << "\t-bed\t"          << "When using BAM input (-abam), write output as BED. The default" << endl;
    cerr                        << "\t\tis to write output in BAM when using -abam." << endl << endl;

    cerr << "\t-s\t"            << "Require same strandedness.  That is, only give records" << endl;
    cerr                        << "\t\tthat have the same strand. Use '-s forward' or '-s reverse'" << endl;
    cerr						<< "\t\tfor forward or reverse strand records, respectively." << endl;
    cerr                        << "\t\t- By default, records are reported without respect to strand." << endl << endl;

    cerr << "\t-header\t"       << "Print the header from the input file prior to results." << endl << endl;

    cerr << "Notes: " << endl;
    cerr << "\tTBD: Enter other usage notes here." << endl << endl;

    // end the program here
    exit(1);

}


