/*****************************************************************************
  bedtools.cpp

  bedtools command line interface.  
  Thanks to Heng Li, as this interface is inspired and 
  based upon his samtools interface.

  (c) 2009-2011 - Aaron Quinlan
  Quinlan Laboratory
  Department of Public Health Sciences
  Center for Public Health genomics
  University of Virginia
  aaronquinlan@gmail.com
  
  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools"

// colors for the term's menu 
#define RESET "\033[m"
#define GREEN "\033[1;32m"
#define BLUE "\033[1;34m"
#define RED "\033[1;31m"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

int annotate_main(int argc, char* argv[]);//
int bamtobed_main(int argc, char* argv[]);//
int bamtofastq_main(int argc, char* argv[]);//
int bed12tobed6_main(int argc, char* argv[]); //
int bedtobam_main(int argc, char* argv[]);//
int bedtoigv_main(int argc, char* argv[]);//
int bedpetobam_main(int argc, char* argv[]);//
int closest_main(int argc, char* argv[]); //
int cluster_main(int argc, char* argv[]); //
int complement_main(int argc, char* argv[]);//
int coverage_main(int argc, char* argv[]); //
int expand_main(int argc, char* argv[]);//
int fastafrombed_main(int argc, char* argv[]);//
int flank_main(int argc, char* argv[]); //
int genomecoverage_main(int argc, char* argv[]);//
int getoverlap_main(int argc, char* argv[]);//
int groupby_main(int argc, char* argv[]);//
int intersect_main(int argc, char* argv[]); //
int links_main(int argc, char* argv[]);//
int maskfastafrombed_main(int argc, char* argv[]);//
int map_main(int argc, char* argv[]); //
int merge_main(int argc, char* argv[]); //
int multibamcov_main(int argc, char* argv[]);//
int multiintersect_main(int argc, char* argv[]);//
int nuc_main(int argc, char* argv[]);//
int pairtobed_main(int argc, char* argv[]);//
int pairtopair_main(int argc, char* argv[]);//
int random_main(int argc, char* argv[]); //
int shuffle_main(int argc, char* argv[]); //
int slop_main(int argc, char* argv[]); //
int sort_main(int argc, char* argv[]); //
int subtract_main(int argc, char* argv[]); //
int tagbam_main(int argc, char* argv[]);//
int unionbedgraphs_main(int argc, char* argv[]);//
int window_main(int argc, char* argv[]); //
int windowmaker_main(int argc, char* argv[]); //
int bedtools_help(void);
int bedtools_faq(void);


int main(int argc, char *argv[])
{
    // make sure the user at least entered a sub_command
    if (argc < 2) return bedtools_help();

    std::string sub_cmd = argv[1];

    // genome arithmetic tools
    if (sub_cmd == "intersect")        return intersect_main(argc-1, argv+1);
    else if (sub_cmd == "window")      return window_main(argc-1, argv+1);
    else if (sub_cmd == "closest")     return closest_main(argc-1, argv+1);
    else if (sub_cmd == "coverage")    return coverage_main(argc-1, argv+1);
    else if (sub_cmd == "map")         return map_main(argc-1, argv+1);
    else if (sub_cmd == "genomecov")   return genomecoverage_main(argc-1, argv+1);
    else if (sub_cmd == "merge")       return merge_main(argc-1, argv+1);
    else if (sub_cmd == "cluster")     return cluster_main(argc-1, argv+1);    
    else if (sub_cmd == "complement")  return complement_main(argc-1, argv+1);
    else if (sub_cmd == "subtract")    return subtract_main(argc-1, argv+1);
    else if (sub_cmd == "slop")        return slop_main(argc-1, argv+1);
    else if (sub_cmd == "flank")       return flank_main(argc-1, argv+1);
    else if (sub_cmd == "sort")        return sort_main(argc-1, argv+1);
    else if (sub_cmd == "random")      return random_main(argc-1, argv+1);
    else if (sub_cmd == "shuffle")     return shuffle_main(argc-1, argv+1);
    else if (sub_cmd == "annotate")    return annotate_main(argc-1, argv+1);

    // Multi-way file comparisonstools
    else if (sub_cmd == "multiinter")  return multiintersect_main(argc-1, argv+1);
    else if (sub_cmd == "unionbedg")   return unionbedgraphs_main(argc-1, argv+1);

    // paired-end conversion tools
    else if (sub_cmd == "pairtobed")   return pairtobed_main(argc-1, argv+1);
    else if (sub_cmd == "pairtopair")  return pairtopair_main(argc-1, argv+1);

    // format conversion tools
    else if (sub_cmd == "bamtobed")    return bamtobed_main(argc-1, argv+1);
    else if (sub_cmd == "bedtobam")    return bedtobam_main(argc-1, argv+1);
    else if (sub_cmd == "bamtofastq")  return bamtofastq_main(argc-1, argv+1);
    else if (sub_cmd == "bedpetobam")  return bedpetobam_main(argc-1, argv+1);
    else if (sub_cmd == "bed12tobed6") return bed12tobed6_main(argc-1, argv+1);

    // BAM-specific tools
    else if (sub_cmd == "multicov")    return multibamcov_main(argc-1, argv+1);
    else if (sub_cmd == "tag")         return tagbam_main(argc-1, argv+1);

    // fasta tools
    else if (sub_cmd == "getfasta")    return fastafrombed_main(argc-1, argv+1);
    else if (sub_cmd == "maskfasta")   return maskfastafrombed_main(argc-1, argv+1);
    else if (sub_cmd == "nuc")         return nuc_main(argc-1, argv+1);

    // misc. tools
    else if (sub_cmd == "overlap")     return getoverlap_main(argc-1, argv+1);
    else if (sub_cmd == "igv")         return bedtoigv_main(argc-1, argv+1);
    else if (sub_cmd == "links")       return links_main(argc-1, argv+1);
    else if (sub_cmd == "makewindows") return windowmaker_main(argc-1, argv+1);
    else if (sub_cmd == "groupby")     return groupby_main(argc-1, argv+1);
    else if (sub_cmd == "expand")      return expand_main(argc-1, argv+1);

    // help
    else if (sub_cmd == "-h" || sub_cmd == "--help" ||
             sub_cmd == "-help")
        return bedtools_help();

    // frequently asked questions
    else if (sub_cmd == "--FAQ" || sub_cmd == "--faq" ||
             sub_cmd == "-FAQ"  || sub_cmd == "-faq")
        return bedtools_faq();

    // verison information
    else if (sub_cmd == "-version" || sub_cmd == "--version")
        cout << "bedtools " << VERSION << endl;

    // verison information
    else if (sub_cmd == "-contact" || sub_cmd == "--contact")
    {
        cout << endl;
        cout << "- For further help, or to report a bug, please " << endl;
        cout << "  email the bedtools mailing list: " << endl;
        cout << "     bedtools-discuss@googlegroups.com" << endl << endl;

        cout << "- Stable releases of bedtools can be found at: " << endl;
        cout << "     http://bedtools.googlecode.com" << endl << endl;

        cout << "- The development repository can be found at: " << endl;
        cout << "     https://github.com/arq5x/bedtools" << endl << endl;
    }
    // unknown
    else {
        // TODO: Implement a Levenstein-based "did you mean???"
        cerr << "error: unrecognized command: " << argv[1] << endl << endl;
        return 1;
    }
    return 0;
}

int bedtools_help(void)
{
    cout  << PROGRAM_NAME  << ": flexible tools for genome arithmetic and DNA sequence analysis.\n";
    cout << "usage:    bedtools <subcommand> [options]" << endl << endl;

    cout  << "The bedtools sub-commands include:" << endl;
    
    cout  << endl;
    cout  << "[ Genome arithmetic ]" << endl;
    cout  << "    intersect     "  << "Find overlapping intervals in various ways.\n";
    cout  << "    window        "  << "Find overlapping intervals within a window around an interval.\n";
    cout  << "    closest       "  << "Find the closest, potentially non-overlapping interval.\n";    
    cout  << "    coverage      "  << "Compute the coverage over defined intervals.\n";
    cout  << "    map           "  << "Apply a function to a column for each overlapping interval.\n";
    cout  << "    genomecov     "  << "Compute the coverage over an entire genome.\n";
    cout  << "    merge         "  << "Combine overlapping/nearby intervals into a single interval.\n";
    cout  << "    cluster       "  << "Cluster (but don't merge) overlapping/nearby intervals.\n";
    cout  << "    complement    "  << "Extract intervals _not_ represented by an interval file.\n";
    cout  << "    subtract      "  << "Remove intervals based on overlaps b/w two files.\n";
    cout  << "    slop          "  << "Adjust the size of intervals.\n";
    cout  << "    flank         "  << "Create new intervals from the flanks of existing intervals.\n";
    cout  << "    sort          "  << "Order the intervals in a file.\n";
    cout  << "    random        "  << "Generate random intervals in a genome.\n";
    cout  << "    shuffle       "  << "Randomly redistrubute intervals in a genome.\n";
    cout  << "    annotate      "  << "Annotate coverage of features from multiple files.\n";
    
    cout  << endl;
    cout  << "[ Multi-way file comparisons ]" << endl;
    cout  << "    multiinter    "  << "Identifies common intervals among multiple interval files.\n";
    cout  << "    unionbedg     "  << "Combines coverage intervals from multiple BEDGRAPH files.\n";

    cout  << endl;
    cout  << "[ Paired-end manipulation ]" << endl;
    cout  << "    pairtobed     "  << "Find pairs that overlap intervals in various ways.\n";
    cout  << "    pairtopair    "  << "Find pairs that overlap other pairs in various ways.\n";

    cout  << endl;
    cout  << "[ Format conversion ]" << endl;
    cout  << "    bamtobed      "  << "Convert BAM alignments to BED (& other) formats.\n";
    cout  << "    bedtobam      "  << "Convert intervals to BAM records.\n";
    cout  << "    bedtofastq    "  << "Convert BAM records to FASTQ records.\n";
    cout  << "    bedpetobam    "  << "Convert BEDPE intervals to BAM records.\n";    
    cout  << "    bed12tobed6   "  << "Breaks BED12 intervals into discrete BED6 intervals.\n";

    cout  << endl;
    cout  << "[ Fasta manipulation ]" << endl;
    cout  << "    getfasta      "  << "Use intervals to extract sequences from a FASTA file.\n";
    cout  << "    maskfasta     "  << "Use intervals to mask sequences from a FASTA file.\n";
    cout  << "    nuc           "  << "Profile the nucleotide content of intervals in a FASTA file.\n";

    cout  << endl;
    cout  << "[ BAM focused tools ]" << endl;
    cout  << "    multicov      "  << "Counts coverage from multiple BAMs at specific intervals.\n";
    cout  << "    tag           "  << "Tag BAM alignments based on overlaps with interval files.\n";

    cout  << endl;
    cout  << "[ Miscellaneous tools ]" << endl;
    cout  << "    overlap       "  << "Computes the amount of overlap from two intervals.\n"; 
    cout  << "    igv           "  << "Create an IGV snapshot batch script.\n";
    cout  << "    links         "  << "Create a HTML page of links to UCSC locations.\n";
    cout  << "    makewindows   "  << "Make interval \"windows\" across a genome.\n";
    cout  << "    groupby       "  << "Group by common cols. & summarize oth. cols. (~ SQL \"groupBy\")\n";
    cout  << "    expand        "  << "Replicate lines based on lists of values in columns.\n";

    cout  << endl;
    cout  << "[ General help ]" << endl;
    cout  << "    --help        "  << "Print this help menu.\n";
    //cout  << "    --faq         "  << "Frequently asked questions.\n";  TODO
    cout  << "    --version     "  << "What version of bedtools are you using?.\n";
    cout  << "    --contact     "  << "Feature requests, bugs, mailing lists, etc.\n";

    cout << "\n";
    return 0;
}


int bedtools_faq(void)
{
    cout << "\n";

    cout << "Q1. How do I see the help for a given command?" << endl;
    cout << "A1. All BEDTools commands have a \"-h\" option. Additionally, some tools " << endl;
    cout << "    will provide the help menu if you just type the command line " << endl;
    cout << "    followed by enter. " << endl;

    cout << "\n";
    return 0;
}
