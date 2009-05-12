#include "windowBed.h"
#include "version.h"

using namespace std;


// define the version
#define PROGRAM_NAME "windowBed"

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
	int leftSlop = 1000;
	int rightSlop = 1000;

	bool haveBedA = false;
	bool haveBedB = false;
	bool noHit = false;
	bool anyHit = false;
	bool writeCount = false;
	bool haveSlop = false;
	bool haveLeft = false;
	bool haveRight = false;
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
		else if(PARAMETER_CHECK("-c", 2, parameterLength)) {
			writeCount = true;
		}
		else if (PARAMETER_CHECK("-v", 2, parameterLength)) {
			noHit = true;
		}
		else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
			forceStrand = true;
		}
		else if (PARAMETER_CHECK("-w", 2, parameterLength)) {
			haveSlop = true;
			leftSlop = atoi(argv[i + 1]);
			rightSlop = leftSlop;
			i++;
		}
		else if (PARAMETER_CHECK("-l", 2, parameterLength)) {
			haveLeft = true;
			leftSlop = atoi(argv[i + 1]);
			i++;
		}
		else if (PARAMETER_CHECK("-r", 2, parameterLength)) {
			haveRight = true;
			rightSlop = atoi(argv[i + 1]);
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
	
	if (anyHit && writeCount) {
		cerr << endl << "*****" << endl << "*****ERROR: Request either -u OR -c, not both." << endl << "*****" << endl;
		showHelp = true;
	}

	if (haveLeft && (leftSlop < 0)) {
		cerr << endl << "*****" << endl << "*****ERROR: Upstream window (-l) must be positive." << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (haveRight && (rightSlop < 0)) {
		cerr << endl << "*****" << endl << "*****ERROR: Downstream window (-r) must be positive." << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (haveSlop && (haveLeft || haveRight)) {
		cerr << endl << "*****" << endl << "*****ERROR: Cannot choose -w with -l or -r.  Either specify -l and -r or specify solely -w" << endl << "*****" << endl;
		showHelp = true;		
	}
	
	if ((haveLeft && !haveRight) || (haveRight && !haveLeft)) {
		cerr << endl << "*****" << endl << "*****ERROR: Please specify both -l and -r." << endl << "*****" << endl;
		showHelp = true;		
	}
	
	if (!showHelp) {
		BedWindow *bi = new BedWindow(bedAFile, bedBFile, leftSlop, rightSlop, anyHit, noHit, writeCount, forceStrand);
		bi->WindowIntersectBed();
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
	cerr << "Description: Examines a \"window\" around each feature in A.bed and" << endl;
	cerr << "reports all features in B.bed that are within the window. For each" << endl;
	cerr << "overlap the entire entry in A and B are reported." << endl << endl;

	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -a <a.bed> -b <b.bed>" << endl << endl;

	cerr << "OPTIONS: " << endl;
	cerr << "\t" << "-s\t\t\t"            	<< "Force strandedness.  Only report hits in B that overlap A on the same strand." << endl << "\t\t\t\tBy default, overlaps are reported without respect to strand." << endl << endl;	
	cerr << "\t" << "-w (def. 1000)\t\t"	<< "Base pairs added upstream and downstream of each entry in A when searching for overlaps in B." << endl << endl;	
	cerr << "\t" << "-l (def. 1000)\t\t"	<< "Base pairs added upstream (left of) of each entry in A when searching for overlaps in B." << endl << endl;	
	cerr << "\t" << "-r (def. 1000)\t\t"	<< "Base pairs added downstream (right of) of each entry in A when searching for overlaps in B." << endl << endl;	
	cerr << "\t" << "-u\t\t\t"            	<< "Write ORIGINAL a.bed entry ONCE if ANY overlap with B.bed." << endl << "\t\t\t\tIn other words, just report the fact >=1 hit was found." << endl << endl;
	cerr << "\t" << "-v \t\t\t"             << "Only report those entries in A that have NO OVERLAP in B within the requested window." << endl << "\t\t\t\tSimilar to grep -v." << endl << endl;
	cerr << "\t" << "-c \t\t\t"				<< "For each entry in A, report the number of hits in B within the requested window." << endl << "\t\t\t\tReports 0 for A entries that have no overlap with B." << endl << endl;

	cerr << "NOTES: " << endl;
	cerr << "\t" << "-i stdin\t\t"	<< "Allows BED file A to be read from stdin.  E.g.: cat a.bed | windowBed -a stdin -b B.bed" << endl << endl;
	cerr << "\t***Only BED3 - BED6 formats allowed.***"<< endl << endl;


	// end the program here
	exit(1);

}
