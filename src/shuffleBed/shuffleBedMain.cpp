#include "shuffleBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "shuffleBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	// input files
	string bedFile;
	string excludeFile;	
	string genomeFile;
	
	bool haveBed = false;
	bool haveGenome = false;
	bool haveExclude = false;
	bool haveSeed = false;
	int seed = -1;
	bool sameChrom = false;

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
		else if(PARAMETER_CHECK("-x", 2, parameterLength)) {
			haveExclude = true;
			excludeFile = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
			haveSeed = true;
			seed = atoi(argv[i + 1]);
			i++;
		}	
		else if(PARAMETER_CHECK("-sameChrom", 10, parameterLength)) {
			sameChrom = true;
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
		BedShuffle *bc = new BedShuffle(bedFile, genomeFile, excludeFile, haveSeed, haveExclude, sameChrom, seed);
		
		if (haveExclude) {
			bc->ShuffleWithExclusions();
		}
		else {
			bc->Shuffle();
		}
		
		return 0;
	}
	else {
		ShowHelp();
	}
}

void ShowHelp(void) {
	
	cerr << "===============================================" << endl;
	cerr << " " <<PROGRAM_NAME << " v" << VERSION << endl ;
	cerr << " Aaron Quinlan, Ph.D. (aaronquinlan@gmail.com)  " << endl ;
	cerr << " Hall Laboratory, University of Virginia" << endl;
	cerr << "===============================================" << endl << endl;
	cerr << "Description: Randomly shuffle the locations of a BED (-i) file among a genome." << endl << endl;
	cerr << "***NOTE: Only tab-delimited BED3 - BED6 formats allowed.***"<< endl;

	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -g <genome> -i <bed>" << endl << endl;
	
	cerr << "OPTIONS: " << endl;
	cerr << "\t" << "-x\t\t\t"            	<< "A BED file of coordinates in which features in -i should not be placed (e.g. gaps.bed)." << endl;
	cerr << "\t" << "-sameChrom\t\t"      << "Keep features in -i on the same chromosome.  By default, the chrom and position are randomly chosen." << endl << endl;
	cerr << "\t" << "-s\t\t\t"     		 	<< "Supply an integer seed for the shuffling.  By default, the seed is chosen automatically." << endl << endl;


	cerr << "NOTES: " << endl;
	cerr << "\tThe genome file should tab delimited and structured as follows: <chr><TAB><size>. For example, Mus musculus:" << endl;
	cerr << "\tchr1\t197195432" << endl;
	cerr << "\tchr2\t181748087" << endl;
	cerr << "\t..." << endl;
	cerr << "\tchrY_random\t58682461" << endl << endl;
	
	cerr << "TIPS:" << endl;
	cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract chromosome sizes. For example, H. sapiens:" << endl << endl;
	cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \"select chrom, size from hg18.chromInfo\"  > hg18.genome" 
		<< endl << endl;
		
	
	// end the program here
	exit(1);
}
