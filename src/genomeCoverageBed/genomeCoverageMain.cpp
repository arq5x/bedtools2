#include <iostream>	
#include "genomeCoverageBed.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "coverageBed"

// define the version
#define VERSION "1.1.0"

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
	int max;
	
	bool haveBed = false;
	bool haveGenome = false;
	bool startSites = false;
	bool eachBase = false;

	// check to see if we should print out some help
	if(argc <= 1) showHelp = true;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if(PARAMETER_CHECK("-h", 2, parameterLength) || 
		PARAMETER_CHECK("--help", 5, parameterLength)) {
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
			i++;
		}
		else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
			startSites = true;
			i++;
		}
		else if(PARAMETER_CHECK("-max", 4, parameterLength)) {
			max = atoi(argv[i + 1]);
			i++;
		}
		else {
		  cout << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have both input files
	if (!haveBed || !haveGenome) {
	  cout << endl << "*****" << endl << "*****ERROR: Need both a BED (-i) and a genome ()-g) file. " << endl << "*****" << endl;
	  showHelp = true;
	}
	
	if (!showHelp) {
		
		BedCoverage *bc = new BedCoverage(bedFile, genomeFile, eachBase, startSites, max);
		bc->CoverageBeds();
		
		return 0;
	}
	else {
		ShowHelp();
	}
}

void ShowHelp(void) {
	
	cout << "=======================================================" << endl;
	cout << PROGRAM_NAME << " v" << VERSION << endl ;
	cout << "Aaron Quinlan, Ph.D." << endl;
	cout << "aaronquinlan@gmail.com" << endl ;
	cout << "Hall Laboratory" << endl;
	cout << "Biochemistry and Molecular Genetics" << endl;
	cout << "University of Virginia" << endl; 
	cout << "=======================================================" << endl << endl;
	cout << "Description: Compute the coverage of a BED (-i) file on genome (-g) file." << endl << endl;
	cout << "***NOTE: Only BED3 - BED6 formats allowed.***"<< endl << endl;

	cout << "Usage: " << PROGRAM_NAME << " [OPTIONS] -g <genome> -i <bed>" << endl << endl;
	
	cout << "OPTIONS: " << endl;
	cout << "\t" << "-d\t\t\t"            	<< "Report the depth at each genome position." << endl << "\t\t\t\tDefault behavior is to report a histogram." << endl << endl;
	cout << "\t" << "-s \t\t\t"             << "Report depth based on start sites." << endl << endl;
	cout << "\t" << "-max\t\t\t"            << "Combine all positions with a depth > max into a single bin in the histogram." << endl << endl;
	cout << "NOTES: " << endl;
	cout << "\tThe genome file should tab delimited and structured as follows: <chr><TAB><size>. For example, Mus musculus:" << endl;
	cout << "\tchr1\t197195432" << endl;
	cout << "\tchr2\t181748087" << endl;
	cout << "\t..." << endl;
	cout << "\tchrY_random\t58682461" << endl << endl;
	cout << "TIPS:" << endl;
	cout << "\tOne can use the UCSC Genome Browser's MySQL database to extract chromosome sizes. For example, H. sapiens:" << endl << endl;
	cout << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \"select chrom, size from hg18.chromInfo\"  > hg18.genome" 
		<< endl << endl;
	// end the program here
	exit(1);
}
