#include "coverageBed.h"
#include "version.h"

using namespace std;

// define the version
#define PROGRAM_NAME "coverageBed"

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

	bool haveBedA = false;
	bool haveBedB = false;

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
			haveBedA = true;
			bedAFile = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
			haveBedB = true;
			bedBFile = argv[i + 1];
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

	if (!showHelp) {

		BedGraph *bg = new BedGraph(bedAFile, bedBFile);
		bg->GraphBed();
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

	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -a <a.bed> -b <b.bed>" << endl << endl;

	cerr << "NOTES: " << endl;
	cerr << "\t***Only BED3 - BED6 formats allowed.***"<< endl << endl;

	//cerr << "OPTIONS: " << endl;
	//cerr << "\t" << "-u\t\t\t"            	<< "Write ORIGINAL a.bed entry ONCE if ANY overlap bed." << endl << "\t\t\t\tIn other words, ignore multiple overlaps." << endl << endl;
	//cerr << "\t" << "-v \t\t\t"             << "Only report those entries in A that have NO OVERLAP in B." << endl << "\t\t\t\tSimilar to grep -v." << endl << endl;
	//cerr << "\t" << "-f (e.g. 0.05)\t\t"	<< "Minimum overlap req'd as fraction of a.bed." << endl << "\t\t\t\tDefault is 1E-9 (effectively 1bp)." << endl << endl;
	//cerr << "\t" << "-c \t\t\t"				<< "For each entry in A, report the number of hits in B while restricting to -f." << endl << "\t\t\t\tReports 0 for A entries that have no overlap with B." << endl << endl;
	//cerr << "\t" << "-wb \t\t\t"			<< "Write the entry in B for each overlap." << endl << "\t\t\t\tUseful for knowing _what_ A overlaps. Restricted by -f." << endl << endl;
	// end the program here
	exit(1);

}
