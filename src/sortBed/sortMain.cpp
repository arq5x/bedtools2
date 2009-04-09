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
		else {
		  cout << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have both input files
	if (!haveBed) {
	  cout << endl << "*****" << endl << "*****ERROR: Need -i BED file. " << endl << "*****" << endl;
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
	
	cout << "=======================================================" << endl;
	cout << PROGRAM_NAME << " v" << VERSION << endl ;
	cout << "Aaron Quinlan, Ph.D." << endl;
	cout << "aaronquinlan@gmail.com" << endl ;
	cout << "Hall Laboratory" << endl;
	cout << "Biochemistry and Molecular Genetics" << endl;
	cout << "University of Virginia" << endl; 
	cout << "=======================================================" << endl << endl;
	cout << "Description: Sorts a BED file by chrom, then by start position." << endl << endl;
	cout << "***NOTE: Only BED3 - BED6 formats allowed.***"<< endl << endl;

	cout << "Usage: " << PROGRAM_NAME << " [OPTIONS] -i <input.bed>" << endl << endl;
	// end the program here
	exit(1);
	
}
