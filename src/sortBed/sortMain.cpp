#include <iostream>	
#include "sortBed.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "sortBed"

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
	bool haveBed = false;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if(PARAMETER_CHECK("-h", 2, parameterLength) || 
		PARAMETER_CHECK("--help", 5, parameterLength)) {
			showHelp = true;
		}
	}

	if(showHelp) ShowHelp();

	// was the file passed in via stdin?
	if(argc == 1) {
		haveBed = true;
	}

	// do some parsing (all of these parameters require 2 strings)
	for(int i = 1; i < argc; i++) {

		int parameterLength = (int)strlen(argv[i]);

		if(argv[i]) {
			haveBed = true;
			bedFile = argv[i];
			i++;
		}
		else {
		  cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have both input files
	if (!haveBed) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need -i BED file. " << endl << "*****" << endl;
	  showHelp = true;
	}
	
	if (!showHelp) {
		BedSort *bm = new BedSort(bedFile);
		bm->SortBed();
		return 0;
	}
	else {
		ShowHelp();
	}
}

void ShowHelp(void) {
	
	cerr << "=======================================================" << endl;
	cerr << PROGRAM_NAME << " v" << VERSION << endl ;
	cerr << "Aaron Quinlan, Ph.D." << endl;
	cerr << "aaronquinlan@gmail.com" << endl ;
	cerr << "Hall Laboratory" << endl;
	cerr << "Biochemistry and Molecular Genetics" << endl;
	cerr << "University of Virginia" << endl; 
	cerr << "=======================================================" << endl << endl;
	cerr << "Description: Sorts a BED file by chrom, then by start position." << endl << endl;
	cerr << "***NOTE: Only BED3 - BED6 formats allowed.***"<< endl << endl;

	cerr << "Usage: " << PROGRAM_NAME << " <input.bed>" << endl << endl;
	// end the program here
	exit(1);
	
}
