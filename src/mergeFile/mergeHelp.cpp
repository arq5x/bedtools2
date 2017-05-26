/*****************************************************************************
  mergeMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "CommonHelp.h"
#include "KeyListOps.h"


void merge_help(void) {
    
    cerr << "\nTool:    bedtools merge (aka mergeBed)" << endl;
    cerr << "Version: " << VERSION << "\n";        
    cerr << "Summary: Merges overlapping BED/GFF/VCF entries into a single interval." << endl << endl;

    cerr << "Usage:   " << "bedtools merge" << " [OPTIONS] -i <bed/gff/vcf>" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-s\t"           << "Force strandedness.  That is, only merge features" << endl;
    cerr                       << "\t\tthat are on the same strand." << endl;
    cerr                       << "\t\t- By default, merging is done without respect to strand." << endl << endl;

    cerr << "\t-S\t"		   << "Force merge for one specific strand only." << endl;
    cerr << "\t\t"             << "Follow with + or - to force merge from only" << endl;
    cerr << "\t\t"			   << "the forward or reverse strand, respectively." << endl;
    cerr << "\t\t"			   << "- By default, merging is done without respect to strand." << endl << endl;

    cerr << "\t-d\t"           << "Maximum distance between features allowed for features" << endl;    cerr                       << "\t\tto be merged." << endl;    cerr                       << "\t\t- Def. 0. That is, overlapping & book-ended features are merged." << endl;
    cerr                       << "\t\t- (INTEGER)" << endl;
    cerr                       << "\t\t- Note: negative values enforce the number of b.p. required for overlap." << endl << endl;


    KeyListOpsHelp();
    allToolsCommonHelp();

    cerr << "Notes: " << endl;
    cerr << "\t(1) The input file (-i) file must be sorted by chrom, then start." << endl << endl;
}
