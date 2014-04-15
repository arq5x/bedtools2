/*****************************************************************************
  mergeMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "mergeFile.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools merge"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void merge_help(void);

int merge_main(int argc, char* argv[]) {

    ContextMerge *context = new ContextMerge();
    if (!context->parseCmdArgs(argc, argv, 1) || context->getShowHelp() || !context->isValidState()) {
        if (!context->getErrorMsg().empty()) {
            cerr << context->getErrorMsg() << endl;
        }
        merge_help();
        delete context;
        return 0;
    }
    MergeFile *mergeFile = new MergeFile(context);

    bool retVal = mergeFile->merge();
    delete mergeFile;
    delete context;
    return retVal ? 0 : 1;
}

void merge_help(void) {
    
    cerr << "\nTool:    bedtools merge (aka mergeBed)" << endl;
    cerr << "Version: " << VERSION << "\n";        
    cerr << "Summary: Merges overlapping BED/GFF/VCF entries into a single interval." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-s\t"                     << "Force strandedness.  That is, only merge features" << endl;
    cerr                                 << "\t\tthat are the same strand." << endl;
    cerr                                 << "\t\t- By default, merging is done without respect to strand." << endl << endl;

    cerr << "\t-n\t"                     << "Report the number of BED entries that were merged." << endl;
    cerr                                 << "\t\t- Note: \"1\" is reported if no merging occurred." << endl << endl;


    cerr << "\t-d\t"                     << "Maximum distance between features allowed for features" << endl;
    cerr                                 << "\t\tto be merged." << endl;
    cerr                                 << "\t\t- Def. 0. That is, overlapping & book-ended features are merged." << endl;
    cerr                                 << "\t\t- (INTEGER)" << endl << endl;

    cerr << "\t-nms\t"                   << "Report the names of the merged features separated by commas." << endl;
    cerr                                 << "\t\tChange delim. with -delim." << endl << endl;
    
    cerr << "\t-scores\t"                << "Report the scores of the merged features. Specify one of " << endl;
    cerr                                 << "\t\tthe following options for reporting scores:" << endl;
    cerr                                 << "\t\t  sum, min, max," << endl;
    cerr                                 << "\t\t  mean, median, mode, antimode," << endl;
    cerr                                 << "\t\t  collapse (i.e., print a semicolon-separated list)," << endl;
    cerr                                 << "\t\t- (INTEGER)" << endl << endl;
    
    cerr << "\t-delim\t"                 << "Specify a custom delimiter for the -nms and -scores concat options" << endl;
    cerr                                 << "\t\t- Example: -delim \"|\"" << endl;
    cerr                                 << "\t\t- Default: \",\"." << endl << endl;
    
    
    cerr << "Notes: " << endl;
    cerr << "\t(1) All output, regardless of input type (e.g., GFF or VCF)" << endl;
    cerr << "\t    will in BED format with zero-based starts" << endl << endl;

    cerr << "\t(2) The input file (-i) file must be sorted by chrom, then start." << endl << endl;

    // end the program here
    exit(1);

}
