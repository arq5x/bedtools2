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
#include "version.h"
#include "BedtoolsDriver.h"

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

bool sub_main(const string &subCmd);
void showHelp(const string &subCmd);
void showErrors(const string &errors);

int annotate_main(int argc, char* argv[]);//
int bamtobed_main(int argc, char* argv[]);//
int bamtofastq_main(int argc, char* argv[]);//
int bed12tobed6_main(int argc, char* argv[]); //
int bedtobam_main(int argc, char* argv[]);//
int bedtoigv_main(int argc, char* argv[]);//
int bedpetobam_main(int argc, char* argv[]);//
void closest_help();
int cluster_main(int argc, char* argv[]); //
void complement_help();
void coverage_help();
int regress_test_main(int argc, char **argv); //
int expand_main(int argc, char* argv[]);//
int fastafrombed_main(int argc, char* argv[]);//
int flank_main(int argc, char* argv[]); //
int genomecoverage_main(int argc, char* argv[]);//
int getoverlap_main(int argc, char* argv[]);//
void groupby_help();
void intersect_help();
void map_help();
void jaccard_help(); //
void fisher_help();
int links_main(int argc, char* argv[]);//
int maskfastafrombed_main(int argc, char* argv[]);//
int map_main(int argc, char* argv[]); //
void merge_help();
int multibamcov_main(int argc, char* argv[]);//
int multiintersect_main(int argc, char* argv[]);//
int nuc_main(int argc, char* argv[]);//
int pairtobed_main(int argc, char* argv[]);//
int pairtopair_main(int argc, char* argv[]);
int random_main(int argc, char* argv[]); //
void summary_help();
int reldist_main(int argc, char* argv[]); //
void sample_help();
int shift_main(int argc, char* argv[]); //
int shuffle_main(int argc, char* argv[]); //
int slop_main(int argc, char* argv[]); //
int split_main(int argc, char* argv[]); //
int sort_main(int argc, char* argv[]); //
void spacing_help();
void subtract_help();
int tagbam_main(int argc, char* argv[]);//
int unionbedgraphs_main(int argc, char* argv[]);//
int window_main(int argc, char* argv[]); //
int windowmaker_main(int argc, char* argv[]); //
int bedtools_help(void);
int bedtools_faq(void);

const char* cram_reference = NULL;

int parse_global_args(int argc, char** argv) {
	for(int i = 1; i < argc - 1; i ++) {
		string this_arg = argv[i];
		if (this_arg == "--cram-ref") {
			cram_reference = argv[i + 1];
			for(int j = i + 2; j < argc; j ++) {
				argv[j - 2] = argv[j];
			}
			i --;
			argc -= 2;
		}
	}
	return argc;
}



int main(int argc, char *argv[])
{
	argc = parse_global_args(argc, argv);
    // make sure the user at least entered a sub_command
    if (argc < 2) return bedtools_help();

    string subCmd(argv[1]);
    BedtoolsDriver btDriver;
    if (btDriver.supports(subCmd)) {

		if (btDriver.subMain(argc, argv)) 
        {
            return 0;
		} 
        else if (btDriver.hadError()) 
        {
			showHelp(subCmd);
            showErrors(btDriver.getErrors());
			return 1;
		}
	}

    // genome arithmetic tools
    else if (subCmd == "window")      return window_main(argc-1, argv+1);
    else if (subCmd == "genomecov")   return genomecoverage_main(argc-1, argv+1);
    else if (subCmd == "cluster")     return cluster_main(argc-1, argv+1);
	  else if (subCmd == "shift")       return shift_main(argc-1, argv+1);
    else if (subCmd == "slop")        return slop_main(argc-1, argv+1);
    else if (subCmd == "split")       return split_main(argc-1, argv+1);
    else if (subCmd == "flank")       return flank_main(argc-1, argv+1);
    else if (subCmd == "sort")        return sort_main(argc-1, argv+1);
    else if (subCmd == "random")      return random_main(argc-1, argv+1);
    else if (subCmd == "shuffle")     return shuffle_main(argc-1, argv+1);
    else if (subCmd == "annotate")    return annotate_main(argc-1, argv+1);

    // Multi-way file comparisonstools
    else if (subCmd == "multiinter")  return multiintersect_main(argc-1, argv+1);
    else if (subCmd == "unionbedg")   return unionbedgraphs_main(argc-1, argv+1);

    // paired-end conversion tools
    else if (subCmd == "pairtobed")   return pairtobed_main(argc-1, argv+1);
    else if (subCmd == "pairtopair")  return pairtopair_main(argc-1, argv+1);

    // format conversion tools
    else if (subCmd == "bamtobed")    return bamtobed_main(argc-1, argv+1);
    else if (subCmd == "bedtobam")    return bedtobam_main(argc-1, argv+1);
    else if (subCmd == "bamtofastq")  return bamtofastq_main(argc-1, argv+1);
    else if (subCmd == "bedpetobam")  return bedpetobam_main(argc-1, argv+1);
    else if (subCmd == "bed12tobed6") return bed12tobed6_main(argc-1, argv+1);

    // BAM-specific tools
    else if (subCmd == "multicov")    return multibamcov_main(argc-1, argv+1);
    else if (subCmd == "tag")         return tagbam_main(argc-1, argv+1);

    // fasta tools
    else if (subCmd == "getfasta")    return fastafrombed_main(argc-1, argv+1);
    else if (subCmd == "maskfasta")   return maskfastafrombed_main(argc-1, argv+1);
    else if (subCmd == "nuc")         return nuc_main(argc-1, argv+1);

    // statistics tools
    else if (subCmd == "reldist")     return reldist_main(argc-1, argv+1);

    // misc. tools
    else if (subCmd == "overlap")     return getoverlap_main(argc-1, argv+1);
    else if (subCmd == "igv")         return bedtoigv_main(argc-1, argv+1);
    else if (subCmd == "links")       return links_main(argc-1, argv+1);
    else if (subCmd == "makewindows") return windowmaker_main(argc-1, argv+1);
    else if (subCmd == "expand")      return expand_main(argc-1, argv+1);
    else if (subCmd == "regresstest")  return regress_test_main(argc, argv); //this command does need all the orig args.
    // help
    else if (subCmd == "-h" || subCmd == "--help" ||
             subCmd == "-help")
        return bedtools_help();

    // frequently asked questions
    else if (subCmd == "--FAQ" || subCmd == "--faq" ||
             subCmd == "-FAQ"  || subCmd == "-faq")
        return bedtools_faq();

    // verison information
    else if (subCmd == "-version" || subCmd == "--version")
        cout << "bedtools " << VERSION << endl;

    // verison information
    else if (subCmd == "-contact" || subCmd == "--contact")
    {
        cout << endl;
        cout << "- For further help, or to report a bug, please " << endl;
        cout << "  email the bedtools mailing list: " << endl;
        cout << "     bedtools-discuss@googlegroups.com" << endl << endl;

        cout << "- The development repository can be found at: " << endl;
        cout << "     https://github.com/arq5x/bedtools" << endl << endl;

        cout << "- Stable releases of bedtools can be found at: " << endl;
        cout << "     https://github.com/arq5x/bedtools2/releases" << endl << endl;
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
    cout  << PROGRAM_NAME  << " is a powerful toolset for genome arithmetic." << endl << endl;
    cout << "Version:   " << VERSION << endl;
    cout << "About:     developed in the quinlanlab.org and by many contributors worldwide." << endl;
    cout << "Docs:      http://bedtools.readthedocs.io/" << endl;
    cout << "Code:      https://github.com/arq5x/bedtools2" << endl;
    cout << "Mail:      https://groups.google.com/forum/#!forum/bedtools-discuss" << endl << endl;
    cout << "Usage:     bedtools <subcommand> [options]" << endl << endl;

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
    cout  << "    shift         "  << "Adjust the position of intervals.\n";
    cout  << "    subtract      "  << "Remove intervals based on overlaps b/w two files.\n";
    cout  << "    slop          "  << "Adjust the size of intervals.\n";
    cout  << "    flank         "  << "Create new intervals from the flanks of existing intervals.\n";
    cout  << "    sort          "  << "Order the intervals in a file.\n";
    cout  << "    random        "  << "Generate random intervals in a genome.\n";
    cout  << "    shuffle       "  << "Randomly redistribute intervals in a genome.\n";
    cout  << "    sample        "  << "Sample random records from file using reservoir sampling.\n";   
    cout  << "    spacing       "  << "Report the gap lengths between intervals in a file.\n";     
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
    cout  << "    bamtofastq    "  << "Convert BAM records to FASTQ records.\n";
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
    cout  << "[ Statistical relationships ]" << endl;
    cout  << "    jaccard       "  << "Calculate the Jaccard statistic b/w two sets of intervals.\n";
    cout  << "    reldist       "  << "Calculate the distribution of relative distances b/w two files.\n";
    cout  << "    fisher        "  << "Calculate Fisher statistic b/w two feature files.\n";

    cout  << endl;
    cout  << "[ Miscellaneous tools ]" << endl;
    cout  << "    overlap       "  << "Computes the amount of overlap from two intervals.\n"; 
    cout  << "    igv           "  << "Create an IGV snapshot batch script.\n";
    cout  << "    links         "  << "Create a HTML page of links to UCSC locations.\n";
    cout  << "    makewindows   "  << "Make interval \"windows\" across a genome.\n";
    cout  << "    groupby       "  << "Group by common cols. & summarize oth. cols. (~ SQL \"groupBy\")\n";
    cout  << "    expand        "  << "Replicate lines based on lists of values in columns.\n";
    cout  << "    split         "  << "Split a file into multiple files with equal records or base pairs.\n"; 
    cout  << "    summary       "  << "Statistical summary of intervals in a file.\n"; 

	cout << endl;
	cout << "[ General Parameters ]" << endl;
	cout << "     --cram-ref    " << "Reference used by a CRAM input" << endl;

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

void showHelp(const string &subCmd) {
	if (subCmd == "intersect") {
		intersect_help();
	} else if (subCmd == "map") {
		map_help();
	} else if (subCmd == "closest") {
		closest_help();
	} else if (subCmd == "merge") {
		merge_help();
	} else if (subCmd == "jaccard") {
		jaccard_help();
	} else if (subCmd == "subtract") {
		subtract_help();
	} else if (subCmd == "sample") {
		sample_help();
	} else if (subCmd == "spacing") {
		spacing_help();
	} else if (subCmd == "fisher") {
		fisher_help();
	} else if (subCmd == "coverage") {
		coverage_help();
	} else if (subCmd == "complement") {
		complement_help();
	} else if (subCmd == "groupby") {
		groupby_help();
	} else if (subCmd == "summary") {
        summary_help();
  }
}

void showErrors(const string &errors) 
{
    cerr << endl << endl << errors << endl;
}
