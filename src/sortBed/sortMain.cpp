#include "sortBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "sortBed"


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
	int sortChoices = 0;

	bool sortBySizeAsc = false;
	bool sortBySizeDesc = false;
	bool sortByChromThenSizeAsc = false;
	bool sortByChromThenSizeDesc = false;
	bool sortByChromThenScoreAsc = false;
	bool sortByChromThenScoreDesc = false;
	
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
		else if(PARAMETER_CHECK("-sizeA", 6, parameterLength)) {
			sortBySizeAsc = true;
			sortChoices++;
			i++;
		}
		else if(PARAMETER_CHECK("-sizeD", 6, parameterLength)) {
			sortBySizeDesc = true;
			sortChoices++;
			i++;
		}
		else if(PARAMETER_CHECK("-chrThenSizeA", 13, parameterLength)) {
			sortByChromThenSizeAsc = true;
			sortChoices++;
			i++;
		}
		else if(PARAMETER_CHECK("-chrThenSizeD", 13, parameterLength)) {
			sortByChromThenSizeDesc = true;
			sortChoices++;
			i++;
		}	
		else if(PARAMETER_CHECK("-chrThenScoreA", 14, parameterLength)) {
			sortByChromThenScoreAsc = true;
			sortChoices++;
			i++;
		}
		else if(PARAMETER_CHECK("-chrThenScoreD", 14, parameterLength)) {
			sortByChromThenScoreDesc = true;
			sortChoices++;
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
	if (sortChoices > 1) {
		cerr << endl << "*****" << endl << "*****ERROR: Sorting options are mutually exclusive.  Please choose just one. " << endl << "*****" << endl;
		showHelp = true;
	}


	if (!showHelp) {
		BedSort *bm = new BedSort(bedFile);

		if (sortBySizeAsc) {
			bm->SortBedBySizeAsc();
		}
		else if (sortBySizeDesc) {
			bm->SortBedBySizeDesc();
		}
		else if (sortByChromThenSizeAsc) {
			bm->SortBedByChromThenSizeAsc();
		}
		else if (sortByChromThenSizeDesc) {
			bm->SortBedByChromThenSizeDesc();
		}
		else if (sortByChromThenScoreAsc) {
			bm->SortBedByChromThenScoreAsc();
		}
		else if (sortByChromThenScoreDesc) {
			bm->SortBedByChromThenScoreDesc();
		}
		else {
			bm->SortBed();			
		}
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
	cerr << "Description: Sorts a BED file in various and useful ways." << endl << endl;
	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -i <input.bed>" << endl << endl;

	cerr << "OPTIONS: " << endl;
	cerr << "\t" << "-sizeA\t\t"	<< "Sort the BED file by feature size in ascending order.  Sorts across all chromosomes." << endl << endl;
	cerr << "\t" << "-sizeD\t\t"	<< "Sort the BED file by feature size in descending order.  Sorts across all chromosomes." << endl << endl;
	cerr << "\t" << "-chrThenSizeA\t"	<< "Sort the BED file by chrom (ascending), then feature size in ascending order." << endl << endl;
	cerr << "\t" << "-chrThenSizeD\t"	<< "Sort the BED file by chrom (ascending), then feature size in descending order." << endl << endl;
	cerr << "\t" << "-chrThenScoreA\t"	<< "Sort the BED file by chrom (ascending), then score in ascending order." << endl << endl;
	cerr << "\t" << "-chrThenScoreD\t"	<< "Sort the BED file by chrom (ascending), then scor size in descending order." << endl << endl;
	
	cerr << "NOTES: " << endl;
	cerr << "\t" << "-i stdin\t\t"	<< "Allows BED file A to be read from stdin.  E.g.: cat a.bed | sortBed -i stdin" << endl << endl;
	cerr << "\t***Only tab-delimited BED3 - BED6 formats allowed.***"<< endl << endl;
	// end the program here
	exit(1);

}
