/*****************************************************************************
  closestMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "closestBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "closestBed"

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
	string tieMode = "all";
	
	bool haveBedA = false;
	bool haveBedB = false;
	bool haveTieMode = false;
	bool forceStrand = false;


	// check to see if we should print out some help
	if(argc <= 1) showHelp = true;

	for(int i = 1; i < argc; i++) {
		int parameterLength = (int)strlen(argv[i]);

		if( (PARAMETER_CHECK("-h", 2, parameterLength)) || 
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
		else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
			forceStrand = true;
		}
		else if (PARAMETER_CHECK("-t", 2, parameterLength)) {
			haveTieMode = true;
			tieMode = argv[i + 1];
			i++;
		}	
	}

	// make sure we have both input files
	if (!haveBedA || !haveBedB) {
		cerr << endl << "*****" << endl << "*****ERROR: Need -a and -b files. " << endl << "*****" << endl;
		showHelp = true;
	}
	
	if (haveTieMode && (tieMode != "all") && (tieMode != "first") 
					&& (tieMode != "last")) {
		cerr << endl << "*****" << endl << "*****ERROR: Request \"all\" or \"first\" or \"last\" for Tie Mode (-t)" << endl << "*****" << endl;
		showHelp = true;		
	}
	
	if (!showHelp) {
		BedClosest *bc = new BedClosest(bedAFile, bedBFile, forceStrand, tieMode);
		bc->DetermineBedInput();
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
	cerr << "SUMMARY: For each feature in BED A, finds the closest " << endl;
	cerr << "\t feature (upstream or downstream) in BED B." << endl << endl;

	cerr << "USAGE:   " << PROGRAM_NAME << " [OPTIONS] -a <a.bed> -b <b.bed>" << endl << endl;

	cerr << "OPTIONS: " << endl;
	cerr << "  " << "-s\t"      << "Force strandedness.  That is, find the closest feature in B that" << endl;
	cerr						<< "\toverlaps A on the same strand." << endl;
	cerr						<< "\t- By default, overlaps are reported without respect to strand." << endl << endl;

	cerr << "  " << "-t\t"     	<< "How ties for closest feature are handled.  This occurs when two" << endl;
	cerr 						<< "\tfeatures in B have exactly the same overlap with a feature in A." << endl;
	cerr						<< "\tBy default, all such features in B are reported." << endl;
	cerr						<< "\tHere are all the options:" << endl;
	cerr 						<< "\t- \"all\"  Report all ties (default)." << endl;
	cerr 						<< "\t- \"first\"  Report the first tie that occurred in the B file." << endl;
	cerr 						<< "\t- \"last\"  Report the last tie that occurred in the B file." << endl << endl;

	
	cerr << "NOTES: " << endl;
	cerr << "  Reports \"none\" for chrom and \"-1\" for all other fields when a feature" << endl;
	cerr << "  is not found in B on the same chromosome as the feature in A." << endl;
	cerr << "  E.g. none\t-1\t-1" << endl << endl;
	
	// end the program here
	exit(1);

}
