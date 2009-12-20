#include "intersectBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "intersectBed"


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

	bool haveBedA = false;
	bool haveBedB = false;
	bool noHit = false;
	bool anyHit = false;
	bool writeA = false;	
	bool writeB = false;
	bool writeCount = false;
	bool haveFraction = false;
	bool reciprocalFraction = false;
	bool forceStrand = false;

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
		else if(PARAMETER_CHECK("-r", 2, parameterLength)) {
			reciprocalFraction = true;
		}
		else if (PARAMETER_CHECK("-v", 2, parameterLength)) {
			noHit = true;
		}
		else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
			forceStrand = true;
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

	if (reciprocalFraction && !haveFraction) {
		cerr << endl << "*****" << endl << "*****ERROR: If using -r, you need to define -f." << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (anyHit && writeCount) {
		cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -c, not both." << endl << "*****" << endl;
		showHelp = true;
	}

	if (!showHelp) {

		BedIntersect *bi = new BedIntersect(bedAFile, bedBFile, anyHit, writeA, writeB, overlapFraction, noHit, writeCount, forceStrand, reciprocalFraction);
		bi->DetermineBedInput();
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
	cerr << "Description: Report overlaps between a.bed and b.bed." << endl << endl;

	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -a <a.bed> -b <b.bed>" << endl << endl;

	cerr << "OPTIONS: " << endl;
	cerr << "\t" << "-s\t\t\t"            	<< "Force strandedness.  Only report hits in B that overlap A on the same strand." << endl << "\t\t\t\tBy default, overlaps are reported without respect to strand." << endl << endl;
	cerr << "\t" << "-u\t\t\t"            	<< "Write ORIGINAL a.bed entry ONCE if ANY overlap with B.bed." << endl << "\t\t\t\tIn other words, just report the fact >=1 hit was found." << endl << endl;
	cerr << "\t" << "-v \t\t\t"             << "Only report those entries in A that have NO OVERLAP in B." << endl << "\t\t\t\tSimilar to grep -v." << endl << endl;
	cerr << "\t" << "-f (e.g. 0.05)\t\t"	<< "Minimum overlap req'd as fraction of a.bed." << endl << "\t\t\t\tDefault is 1E-9 (effectively 1bp)." << endl << endl;
	cerr << "\t" << "-r \t\t\t"				<< "Require that the fraction overlap be reciprocal for A and B."   << endl << "\t\t\t\tIn other words, if -f is 0.90 and -r is used, this requires that" << endl << "\t\t\t\tB overlap 90% of A and A _also_ overlaps 90% of B." << endl << endl;
	cerr << "\t" << "-c \t\t\t"				<< "For each entry in A, report the number of hits in B while restricting to -f." << endl << "\t\t\t\tReports 0 for A entries that have no overlap with B." << endl << endl;
	cerr << "\t" << "-wa \t\t\t"			<< "Write the original entry in A for each overlap." << endl << endl;
	cerr << "\t" << "-wb \t\t\t"			<< "Write the original entry in B for each overlap." << endl << "\t\t\t\tUseful for knowing _what_ A overlaps. Restricted by -f." << endl << endl;

	cerr << "NOTES: " << endl;
	cerr << "\t" << "-i stdin\t\t"	<< "Allows BED file A to be read from stdin.  E.g.: cat a.bed | intersectBed -a stdin -b B.bed" << endl << endl;
	cerr << "\t***Only tab-delimited BED3 - BED6 formats allowed.***"<< endl << endl;

	// end the program here
	exit(1);

}
