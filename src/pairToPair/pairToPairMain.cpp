#include "pairToPair.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "peIntersectBed"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	// input files
	string bedAFile;
	string bedBFile;
	
	// input arguments
	float overlapFraction = 1E-9;
	string searchType = "both";

	// flags to track parameters
	bool haveBedA = false;
	bool haveBedB = false;
	bool haveSearchType = false;
	bool haveFraction = false;

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

		if(PARAMETER_CHECK("-a", 2, parameterLength)) {
			if ((i+1) < argc) {
				haveBedA = true;
				bedAFile = argv[i + 1];
			}
			i++;
		}
		else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
			if ((i+1) < argc) {
				haveBedB = true;
				bedBFile = argv[i + 1];
			}
			i++;
		}	
		else if(PARAMETER_CHECK("-type", 5, parameterLength)) {
			if ((i+1) < argc) {
				haveSearchType = true;
				searchType = argv[i + 1];
			}
			i++;
		}
		else if(PARAMETER_CHECK("-f", 2, parameterLength)) {
			haveFraction = true;
			overlapFraction = atof(argv[i + 1]);
			i++;
		}
		else {
			cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}


	// make sure we have both input files
	if (!haveBedA || !haveBedB) {
		cerr << endl << "*****" << endl << "*****ERROR: Need -a and -b files. " << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (haveSearchType && (searchType != "neither") && (searchType != "both")) {
		cerr << endl << "*****" << endl << "*****ERROR: Request \"both\" or \"neither\"" << endl << "*****" << endl;
		showHelp = true;		
	}

	if (!showHelp) {

		PairToPair *bi = new PairToPair(bedAFile, bedBFile, overlapFraction, searchType);
		bi->IntersectPairs();
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
	cerr << "Description: Report overlaps between a.bedpe and b.bedpe." << endl << endl;

	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -a <a.bedpe> -b <b.bedpe>" << endl << endl;

	cerr << "OPTIONS: " << endl;
	cerr << "\t" << "-f\t\t"	    			<< "Minimum overlap req'd as fraction of a.bed (e.g. 0.05)." << endl << "\t\t\tDefault is 1E-9 (effectively 1bp)." << endl << endl;
	cerr << "\t" << "-type \t\t"				<< "neither\t\t\tReport overlaps if _neither_ end of BEDPE (A) overlaps BEDPE B." << endl;
	cerr 	 									<< "\t\t\tboth\t\t\tReport overlaps if _both_ ends of BEDPE (A) overlap BEDPE B." << endl;
	
	
	cerr << "NOTES: " << endl;
	cerr << "\t" << "-i stdin\t"	<< "Allows BEDPE file A to be read from stdin.  E.g.: cat a.bedpe | peIntersectBed -a stdin -b B.bed" << endl << endl;
	cerr << "\t***Only 10 columns BEDPE formats allowed.  See below).***"<< endl << endl;
	
	cerr << "BEDPE FORMAT (tab-delimited): " << endl;
	cerr << "\t" << "1. chrom for end 1" << endl; 
	cerr << "\t" << "2. start for end 1" << endl;
	cerr << "\t" << "3. end for end 1" << endl;
	cerr << "\t" << "4. chrom for end 2" << endl; 
	cerr << "\t" << "5. start for end 2" << endl;
	cerr << "\t" << "6. end for end 2" << endl;
	cerr << "\t" << "7. name" << endl; 
	cerr << "\t" << "8. score" << endl;
	cerr << "\t" << "9. strand for end 1" << endl;
	cerr << "\t" << "10. strand for end 2" << endl;
	cerr << "\t" << "Note: Strands for each end must be provided if you choose to include strand information." << endl << endl;

	cerr << "Example BEDPE record:" << endl;
	cerr << "chr1\t1000000\t1000100\tchr1\t1001000\t1001100\tsomename\t100\t+\t-" << endl;

	// end the program here
	exit(1);

}
