#include "intersectBed.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "intersectBed"

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
	string bedAFile;
	string bedBFile;

	// input arguments
	float overlapFraction = 1E-9;
	int slop = 0;

	bool haveBedA = false;
	bool haveBedB = false;
	bool noHit = false;
	bool anyHit = false;
	bool writeA = false;	
	bool writeB = false;
	bool writeCount = false;
	bool haveFraction = false;
	bool haveSlop = false;

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

		if(PARAMETER_CHECK("-a", 2, parameterLength)) {
			haveBedA = true;
			bedAFile = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
			haveBedB = true;
			bedBFile = argv[i + 1];
			i++;
		}	
		else if(PARAMETER_CHECK("-u", 2, parameterLength)) {
			anyHit = true;
		}
		else if(PARAMETER_CHECK("-f", 2, parameterLength)) {
			haveFraction = true;
			overlapFraction = atof(argv[i + 1]);
			i++;
		}
		else if(PARAMETER_CHECK("-wa", 3, parameterLength)) {
			writeA = true;
		}
		else if(PARAMETER_CHECK("-wb", 3, parameterLength)) {
			writeB = true;
		}
		else if(PARAMETER_CHECK("-c", 2, parameterLength)) {
			writeCount = true;
		}
		else if (PARAMETER_CHECK("-v", 2, parameterLength)) {
			noHit = true;
		}
		else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
			haveSlop = true;
			slop = atoi(argv[i + 1]);
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

	if (anyHit && noHit) {
		cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -v, not both." << endl << "*****" << endl;
		showHelp = true;
	}

	if (writeB && writeCount) {
		cerr << endl << "*****" << endl << "*****ERROR: Request either -wb OR -c, not both." << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (anyHit && writeCount) {
		cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -c, not both." << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (haveSlop && haveFraction) {
		cerr << endl << "*****" << endl << "*****ERROR: Request either -s OR -f, not both." << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (haveSlop && (slop < 0)) {
		cerr << endl << "*****" << endl << "*****ERROR: Slop (-s) must be positive." << endl << "*****" << endl;
		showHelp = true;
	}

	if (!showHelp) {

		BedIntersect *bi = new BedIntersect(bedAFile, bedBFile, anyHit, writeA, writeB, overlapFraction, slop, noHit,  writeCount);
		bi->IntersectBed();
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
	cerr << "Description: Report overlaps between a.bed and b.bed." << endl << endl;
	cerr << "***NOTE: Only BED3 - BED6 formats allowed.***"<< endl << endl;

	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -a <a.bed> -b <b.bed>" << endl << endl;

	cerr << "OPTIONS: " << endl;
	cerr << "\t" << "-u\t\t\t"            	<< "Write ORIGINAL a.bed entry ONCE if ANY overlap bed." << endl << "\t\t\t\tIn other words, ignore multiple overlaps." << endl << endl;
	cerr << "\t" << "-v \t\t\t"             << "Only report those entries in A that have NO OVERLAP in B." << endl << "\t\t\t\tSimilar to grep -v." << endl << endl;
	cerr << "\t" << "-s (100000)\t\t"	    << "Slop added before and after each entry in A" << endl << "\t\t\t\tUseful for finding overlaps within N bases upstream and downstream." << endl << endl;	
	cerr << "\t" << "-f (e.g. 0.05)\t\t"	<< "Minimum overlap req'd as fraction of a.bed." << endl << "\t\t\t\tDefault is 1E-9 (effectively 1bp)." << endl << endl;
	cerr << "\t" << "-c \t\t\t"				<< "For each entry in A, report the number of hits in B while restricting to -f." << endl << "\t\t\t\tReports 0 for A entries that have no overlap with B." << endl << endl;
	cerr << "\t" << "-wa \t\t\t"			<< "Write the original entry in A for each overlap." << endl << endl;
	cerr << "\t" << "-wb \t\t\t"			<< "Write the original entry in B for each overlap." << endl << "\t\t\t\tUseful for knowing _what_ A overlaps. Restricted by -f." << endl << endl;

	cerr << "NOTES: " << endl;
	cerr << "\t" << "-i stdin\t\t"	<< "Allows intersectBed to read BED from stdin.  E.g.: cat a.bed | intersectBed -a stdin -b B.bed" << endl << endl;


	// end the program here
	exit(1);

}
