/*****************************************************************************
  bedtools.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "version.h"

using namespace std;


// define our program name
#define PROGRAM_NAME "bedtools"

// colors for the term's menu 
#define RESET "\e[m"
#define GREEN "\e[1;32m"
#define BLUE "\e[1;34m"
#define RED "\e[1;31m"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

int annotate_main(int argc, char* argv[]);//
int bamtobed_main(int argc, char* argv[]);//
int bed12tobed6_main(int argc, char* argv[]); //
int bedtobam_main(int argc, char* argv[]);//
int bedtoigv_main(int argc, char* argv[]);//
int bedpetobam_main(int argc, char* argv[]);//
int closest_main(int argc, char* argv[]); //
int complement_main(int argc, char* argv[]);//
int coverage_main(int argc, char* argv[]); //
int fastafrombed_main(int argc, char* argv[]);//
int flank_main(int argc, char* argv[]); //
int genomecoverage_main(int argc, char* argv[]);//
int getoverlap_main(int argc, char* argv[]);//
int intersect_main(int argc, char* argv[]); //
int links_main(int argc, char* argv[]);//
int maskfastafrombed_main(int argc, char* argv[]);//
int merge_main(int argc, char* argv[]); //
int multibamcov_main(int argc, char* argv[]);//
int multiintersect_main(int argc, char* argv[]);//
int nuc_main(int argc, char* argv[]);//
int pairtobed_main(int argc, char* argv[]);//
int pairtopair_main(int argc, char* argv[]);//
int shuffle_main(int argc, char* argv[]); //
int slop_main(int argc, char* argv[]); //
int sort_main(int argc, char* argv[]); //
int subtract_main(int argc, char* argv[]); //
int tagbam_main(int argc, char* argv[]);//
int unionbedgraphs_main(int argc, char* argv[]);//
int window_main(int argc, char* argv[]); //
int bedtools_help(void);
int bedtools_faq(void);

/*
  bedtools command line interface.  
  
  Thanks to Heng Li, as this interface is inspired and based upon his samtools interface.
*/

int main(int argc, char *argv[])
{
    // make sure the user at least entered a sub_command
    if (argc < 2) return bedtools_help();

    // genome arithmetic tools
    if (strcmp(argv[1], "intersect") == 0)        return intersect_main(argc-1, argv+1);
    else if (strcmp(argv[1], "window") == 0)      return window_main(argc-1, argv+1);
    else if (strcmp(argv[1], "closest") == 0)     return closest_main(argc-1, argv+1);
    else if (strcmp(argv[1], "coverage") == 0)    return coverage_main(argc-1, argv+1);
    else if (strcmp(argv[1], "genomecov") == 0)   return genomecoverage_main(argc-1, argv+1);
    else if (strcmp(argv[1], "merge") == 0)       return merge_main(argc-1, argv+1);
    else if (strcmp(argv[1], "complement") == 0)  return complement_main(argc-1, argv+1);
    else if (strcmp(argv[1], "subtract") == 0)    return subtract_main(argc-1, argv+1);
    else if (strcmp(argv[1], "slop") == 0)        return slop_main(argc-1, argv+1);
    else if (strcmp(argv[1], "flank") == 0)       return flank_main(argc-1, argv+1);
    else if (strcmp(argv[1], "sort") == 0)        return sort_main(argc-1, argv+1);
    else if (strcmp(argv[1], "shuffle") == 0)     return shuffle_main(argc-1, argv+1);
    else if (strcmp(argv[1], "annotate") == 0)    return annotate_main(argc-1, argv+1);

    // Multi-way file comparisonstools
    else if (strcmp(argv[1], "multiinter") == 0)  return multiintersect_main(argc-1, argv+1);
    else if (strcmp(argv[1], "unionbedg") == 0)   return unionbedgraphs_main(argc-1, argv+1);

    // paired-end conversion tools
    else if (strcmp(argv[1], "pairtobed") == 0)   return pairtobed_main(argc-1, argv+1);
    else if (strcmp(argv[1], "pairtopair") == 0)  return pairtopair_main(argc-1, argv+1);

    // format conversion tools
    else if (strcmp(argv[1], "bamtobed") == 0)    return bamtobed_main(argc-1, argv+1);
    else if (strcmp(argv[1], "bedtobam") == 0)    return bedtobam_main(argc-1, argv+1);
    else if (strcmp(argv[1], "bedpetobam") == 0)  return bedpetobam_main(argc-1, argv+1);
    else if (strcmp(argv[1], "bed12tobed6") == 0) return bed12tobed6_main(argc-1, argv+1);

    // BAM-specific tools
    else if (strcmp(argv[1], "multicov") == 0)    return multibamcov_main(argc-1, argv+1);
    else if (strcmp(argv[1], "tag") == 0)         return tagbam_main(argc-1, argv+1);    

    // fasta tools
    else if (strcmp(argv[1], "getfasta") == 0)       return fastafrombed_main(argc-1, argv+1);
    else if (strcmp(argv[1], "maskfasta") == 0)      return maskfastafrombed_main(argc-1, argv+1);
    else if (strcmp(argv[1], "nuc") == 0)         return nuc_main(argc-1, argv+1);

    // misc. tools
    else if (strcmp(argv[1], "overlap") == 0)     return getoverlap_main(argc-1, argv+1);
    else if (strcmp(argv[1], "igv") == 0)         return bedtoigv_main(argc-1, argv+1);
    else if (strcmp(argv[1], "links") == 0)       return links_main(argc-1, argv+1);    

    // help
    else if (strcmp(argv[1], "-h") == 0)          return bedtools_help();
    else if (strcmp(argv[1], "--help") == 0)      return bedtools_help();

    else if (strcmp(argv[1], "FAQ") == 0)          return bedtools_faq();
    // unknown
    else {
        cerr << "error: unrecognized command: " << argv[1] << endl << endl;
        return 1;
    }
    return 0;   
}

int bedtools_help(void)
{
    cout << "\n";
    cout << RED << PROGRAM_NAME << RESET << ": flexible tools for genome arithmetic and analysis.\n";
    cout << "Version:  " << VERSION << "\n";
    cout << "Authors:  " << "Aaron Quinlan and others (see THANKS)" << "\n\n";
    cout << "Usage:   bedtools <tool> [options]\n";

    cout << "\nGenome arithmetic:" << endl;
    cout << RED << "    intersect     " << RESET << "Find overlapping intervals in various ways.\n";
    cout << RED << "    window        " << RESET << "Find overlapping intervals within a window around an interval.\n";
    cout << RED << "    closest       " << RESET << "Find the closest, potentially non-overlapping interval.\n";    
    cout << RED << "    coverage      " << RESET << "Compute the coverage over defined intervals.\n";
    cout << RED << "    genomecov     " << RESET << "Compute the coverage over an entire genome.\n";
    cout << RED << "    merge         " << RESET << "Combine overlapping/nearby intervals into a single interval.\n";
    cout << RED << "    complement    " << RESET << "Extract intervals _not_ represented by an interval file.\n";
    cout << RED << "    subtract      " << RESET << "Remove intervals based on overlaps b/w two files.\n";
    cout << RED << "    slop          " << RESET << "Adjust the size of intervals.\n";
    cout << RED << "    flank         " << RESET << "Create new intervals from the flanks of existing intervals.\n";    
    cout << RED << "    sort          " << RESET << "Order the intervals in a file.\n";
    cout << RED << "    shuffle       " << RESET << "Randomly redistrubute intervals in a genome.\n";
    cout << RED << "    annotate      " << RESET << "Annotate coverage of features from multiple files.\n";
    
    cout << "\nMulti-way file comparisons:" << endl;
    cout << RED << "    multiinter    " << RESET << "Identifies common intervals among multiple interval files.\n";
    cout << RED << "    unionbedg     " << RESET << "Combines coverage intervals from multiple BEDGRAPH files.\n";

    cout << "\nPaired-end manipulation:" << endl;
    cout << RED << "    pairtobed     " << RESET << "Find pairs that overlap intervals in various ways.\n";
    cout << RED << "    pairtopair    " << RESET << "Find pairs that overlap other pairs in various ways.\n";

    cout << "\nFormat conversion:\n";   
    cout << RED << "    bamtobed      " << RESET << "Convert BAM alignments to BED (& other) formats.\n";
    cout << RED << "    bedtobam      " << RESET << "Convert intervals to BAM records.\n";
    cout << RED << "    bedpetobam    " << RESET << "Convert BEDPE intervals to BAM records.\n";    
    cout << RED << "    bed12tobed6   " << RESET << "Breaks BED12 intervals into discrete BED6 intervals.\n";

    cout << "\nFasta manipulation:\n";  
    cout << RED << "    getfasta      " << RESET << "Use intervals to extract sequences from a FASTA file.\n";
    cout << RED << "    maskfasta     " << RESET << "Use intervals to mask sequences from a FASTA file.\n";
    cout << RED << "    nuc           " << RESET << "Profile the nucleotide content of intervals in a FASTA file.\n";

    cout << "\nBAM focused tools:\n";   
    cout << RED << "    multicov      " << RESET << "Counts coverage from multiple BAMs at specific intervals.\n";
    cout << RED << "    tag           " << RESET << "Tag BAM alignments based on overlaps with interval files.\n";

    cout << "\nMiscellaneous tools:\n"; 
    cout << RED << "    overlap       " << RESET << "Computes the amount of overlap from two intervals.\n"; 
    cout << RED << "    igv           " << RESET << "Create an IGV snapshot batch script.\n";
    cout << RED << "    links         " << RESET << "Create a HTML page of links to UCSC locations.\n";
        
    cout << "\nGeneral help:\n";    
    cout << RED << "    FAQ           " << RESET << "Frequently asked questions.\n";

    cout << "\n";
    return 1;
}


int bedtools_faq(void)
{
    cout << "\n";
    
    cout << "Q1. How do I see the help for a given command?" << endl;
    cout << "A1. All BEDTools commands have a \"-h\" option. Additionally, some tools " << endl;
    cout << "    will provide the help menu if you just type the command line " << endl;
    cout << "    followed by enter. " << endl;
    
    cout << "Q2. How do I see the help for a given command?" << endl;
    cout << "A1. All BEDTools commands have a \"-h\" option. Additionally, some tools " << endl;
    cout << "    will provide the help menu if you just type the command line " << endl;
    cout << "    followed by enter. " << endl;
    
    
    cout << "\n";
    return 1;
}