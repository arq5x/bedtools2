#include "maskFastaFromBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "maskFastaFromBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	// input files
	string fastaInFile;
	string bedFile;

	// output files
	string fastaOutFile;

	// checks for existence of parameters
	bool haveFastaIn = false;
	bool haveBed = false;
	bool haveFastaOut = false;
	bool softMask = false;

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

		if(PARAMETER_CHECK("-in", 3, parameterLength)) {
			haveFastaIn = true;
			fastaInFile = argv[i + 1];
			i++;
		} 
		else if(PARAMETER_CHECK("-fo", 3, parameterLength)) {
			haveFastaOut = true;
			fastaOutFile = argv[i + 1];
			i++;
		} 
		else if(PARAMETER_CHECK("-bed", 4, parameterLength)) {
			haveBed = true;
			bedFile = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-soft", 5, parameterLength)) {
			softMask = true;
			i++;
		}	
		else {
			cerr << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	if (!haveFastaIn || !haveFastaOut || !haveBed) {
		showHelp = true;
	}
	
	if (!showHelp) {

		MaskFastaFromBed *maskFasta = new MaskFastaFromBed(fastaInFile, bedFile, fastaOutFile, softMask);
		maskFasta->MaskFasta(); 
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
	cerr << "Description: Mask a fasta file based on BED coordinates." << endl << endl;

	cerr << "Usage: " << PROGRAM_NAME << " -in -out -bed (-soft)" << endl << endl;

	cerr << "OPTIONS:" << endl;
	cerr << "\t-in\t\t<FASTA input file>          " << endl;
	cerr << "\t-bed\t\t<BED coordinates file>          " << endl;
	cerr << "\t-fo\t\t<FASTA output file>                           " << endl;
	cerr << "\t-soft\t\tEnforce \"soft\" masking.  That is, instead of masking with Ns," << endl;
	cerr << "\t\t\tmask with lower-case bases." << endl;
	
	cerr << "NOTES: " << endl;
	cerr << "\t***Only tab-delimited BED3 - BED6 formats allowed.***"<< endl << endl;

	// end the program here
	exit(1);
	
}
