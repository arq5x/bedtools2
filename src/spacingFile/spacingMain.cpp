/*
 * spacingMain.cpp
 *
 *  Created on: Jan 15, 2015
 *      Author: Aaron Quinlan
 */

#include <iostream>
#include "ContextSpacing.h"
#include "SpacingFile.h"

#define PROGRAM_NAME "bedtools spacing"

void spacing_help(void);

int spacing_main(int argc, char **argv)
{
    ContextSpacing *context = new ContextSpacing();
    if (!context->parseCmdArgs(argc, argv, 1) || context->getShowHelp() || !context->isValidState()) {
        if (!context->getErrorMsg().empty()) {
            cerr << context->getErrorMsg() << endl;
        }
        spacing_help();
        delete context;
        return 1;
    }
    SpacingFile *spacingFile = new SpacingFile(context);

    bool retVal = spacingFile->getSpacing();
    delete spacingFile;
    delete context;
    return retVal ? 0 : 1;

}

void spacing_help(void) {

    cerr << "\nTool:    bedtools spacing" << endl;
    cerr << "Version: " << VERSION << "\n";
    cerr << "Summary: Report (last col.) the gap lengths between intervals in a file." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf/bam>" << endl << endl;
    
    cerr << "Notes: " << endl;
    cerr << "\t(1)  Input must be sorted by chrom,start (sort -k1,1 -k2,2n for BED)." << endl;
    cerr << "\t(2)  The 1st element for each chrom will have NULL distance. (\".\")." << endl;
    cerr << "\t(3)  Distance for overlapping intervaks is -1 and bookended is 0." << endl << endl;

    cerr << "Example: " << endl;
    cerr << "\t$ cat test.bed " << endl;
    cerr << "\tchr1    0   10 " << endl;
    cerr << "\tchr1    10  20 " << endl;
    cerr << "\tchr1    21  30 " << endl;
    cerr << "\tchr1    35  45 " << endl;
    cerr << "\tchr1    100 200 " << endl << endl;

    cerr << "\t$ bedtools spacing -i test.bed " << endl;
    cerr << "\tchr1    0   10  . " << endl;
    cerr << "\tchr1    10  20  0 " << endl;
    cerr << "\tchr1    21  30  1 " << endl;
    cerr << "\tchr1    35  45  5 " << endl;
    cerr << "\tchr1    100 200 55 " << endl << endl;

    // end the program here
    exit(1);
}
