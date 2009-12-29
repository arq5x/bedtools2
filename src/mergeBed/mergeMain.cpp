/*****************************************************************************
  mergeMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "mergeBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "mergeBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	// input files
	string bedFile;
	int maxDistance = 0;

	// input arguments
	bool numEntries = false;
	bool haveMaxDistance = false;
	bool haveBed = false;
	bool forceStrand = false;
	bool reportNames = false;

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
		else if(PARAMETER_CHECK("-n", 2, parameterLength)) {
			numEntries = true;
		}
		else if(PARAMETER_CHECK("-d", 2, parameterLength)) {
			haveMaxDistance = true;
			maxDistance = atoi(argv[i + 1]);
			i++;
		}
		else if (PARAMETER_CHECK("-s", 2, parameterLength)) {
			forceStrand = true;
		}
		else if (PARAMETER_CHECK("-nms", 4, parameterLength)) {
			reportNames = true;
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
		BedMerge *bm = new BedMerge(bedFile, numEntries, maxDistance, forceStrand, reportNames);
		
		if (!forceStrand) {
			bm->MergeBed();
		}
		else {
			bm->MergeBedStranded();			
		}
		return 0;
	}
	else {
		ShowHelp();
	}
}

void ShowHelp(void) {

	cerr << endl << "PROGRAM: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl << endl;
	
	cerr << "AUTHOR:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl << endl ;
	
	cerr << "SUMMARY: Merges overlapping BED entries into a single interval." << endl << endl;

	cerr << "USAGE:   " << PROGRAM_NAME << " [OPTIONS] -i <input.bed>" << endl << endl;

	cerr << "OPTIONS: " << endl;
	cerr << "  " << "-s\t"      << "Force strandedness.  That is, only merge features" << endl;
	cerr						<< "\tthat are the same strand." << endl;
	cerr						<< "\t- By default, merging is done without respect to strand." << endl << endl;

	cerr << "  " << "-n\t"		<< "Report the number of BED entries that were merged." << endl;
	cerr						<< "\t- Note: \"1\" is reported if no merging occured." << endl << endl;


	cerr << "  " << "-d\t"  	<< "Maximum distance between features allowed for features to be merged." << endl;
	cerr 	 					<< "\t- Def. 0. That is, overlapping and/or book-ended features are merged." << endl;
	cerr						<< "\t- INTEGER" << endl << endl;
	
	cerr << "  " << "-nms\t"  	<< "Report the names of the merged features separated by semicolons." << endl << endl;
	

	// end the program here
	exit(1);

}
