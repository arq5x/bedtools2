/*****************************************************************************
  genomeCoverageMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "genomeCoverageBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "genomeCoverageBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	// input files
	string bedFile;
	string genomeFile;
	int max = 999999999;
	
	bool haveBed = false;
	bool haveGenome = false;
	bool startSites = false;
	bool eachBase = false;

	// check to see if we should print out some help
	if(argc <= 1) showHelp = true;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if((PARAMETER_CHECK("-h", 2, parameterLength)) || 
		(PARAMETER_CHECK("--help", 5, parameterLength))) {
			showHelp = true;
		}
	}

	if(showHelp) ShowHelp();

	// do some parsing (all of these parameters require 2 strings)
	for(int i = 1; i < argc; i++) {

		int parameterLength = (int)strlen(argv[i]);

 		if(PARAMETER_CHECK("-i", 2, parameterLength)) {
			haveBed = true;
			bedFile = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-g", 2, parameterLength)) {
			haveGenome = true;
			genomeFile = argv[i + 1];
			i++;
		}	
		else if(PARAMETER_CHECK("-d", 2, parameterLength)) {
			eachBase = true;
		}
		//else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
		//	startSites = true;
		//	i++;
		//}
		else if(PARAMETER_CHECK("-max", 4, parameterLength)) {
			max = atoi(argv[i + 1]);
			i++;
		}
		else {
		  cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have both input files
	if (!haveBed || !haveGenome) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need both a BED (-i) and a genome (-g) file. " << endl << "*****" << endl;
	  showHelp = true;
	}
	
	if (!showHelp) {
		
		BedCoverage *bc = new BedCoverage(bedFile, genomeFile, eachBase, startSites, max);
		bc->DetermineBedInput();
		
		return 0;
	}
	else {
		ShowHelp();
	}
}

void ShowHelp(void) {
	
	cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

	cerr << "Summary: Compute the coverage of a BED file among a genome." << endl << endl;

	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -g <genome> -i <bed>" << endl << endl;
	
	cerr << "Options: " << endl;
	cerr << "\t-d\t\t"	     	<< "Report the depth at each genome position." << endl;
	cerr 						<< "\t\t\tDefault behavior is to report a histogram." << endl << endl;

	cerr << "\t-max\t"          << "Combine all positions with a depth >= max into" << endl;
	cerr						<< "\t\t\ta single bin in the histogram." << endl;
	cerr						<< "\t\t\t- (INTEGER)" << endl << endl;

	cerr << "Notes: " << endl;
	cerr << "\t(1)  The genome file should tab delimited and structured as follows:" << endl;
	cerr << "\t     <chromName><TAB><chromSize>" << endl << endl;
	cerr << "\tFor example, Human (hg19):" << endl;
	cerr << "\tchr1\t249250621" << endl;
	cerr << "\tchr2\t243199373" << endl;
	cerr << "\t..." << endl;
	cerr << "\tchr18_gl000207_random\t4262" << endl << endl;

	cerr << "\t(2)  NOTE: The input BED file must be grouped by chromosome." << endl;
	cerr << "\t     A simple \"sort -k 1,1 <BED> > <BED>.sorted\" will suffice."<< endl << endl;

	
	cerr << "Tips: " << endl;
	cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract" << endl;
	cerr << "\tchromosome sizes. For example, H. sapiens:" << endl << endl;
	cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e /" << endl;
	cerr << "\t\"select chrom, size from hg19.chromInfo\"  > hg19.genome" << endl << endl;
		
	
	// end the program here
	exit(1);
}
