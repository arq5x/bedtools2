/*****************************************************************************
  pairToBedMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "pairToBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "pairToBed"

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
	string searchType = "either";

	// flags to track parameters
	bool haveBedA = false;
	bool haveBedB = false;
	bool haveSearchType = false;
	bool haveFraction = false;
	bool forceStrand = false;
	bool inputIsBam = false;
	bool outputIsBam = false;
	
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
		else if(PARAMETER_CHECK("-abam", 5, parameterLength)) {
			if ((i+1) < argc) {
				haveBedA = true;
				inputIsBam = true;
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
		else if(PARAMETER_CHECK("-outbam", 7, parameterLength)) {
			outputIsBam = true;
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
	
	if (haveSearchType && (searchType != "either") && (searchType != "neither") && (searchType != "both") 
					    && (searchType != "xor") && (searchType != "ispan") && (searchType != "ospan")
						&& (searchType != "notispan") && (searchType != "notospan")) {
		cerr << endl << "*****" << endl << "*****ERROR: Request \"either\" or \"both\" or \"neither\" or \"xor\" or \"ispan\" or \"ospan\" or \"notispan\" or \"notospan\"" << endl << "*****" << endl;
		showHelp = true;		
	}
	
	if ( ((searchType == "ispan") || (searchType == "ospan") || (searchType == "notispan") || (searchType == "notospan")) 
	     && forceStrand ) {
		cerr << endl << "*****" << endl << "*****ERROR: Cannot enforce strandedness with selected searchtype" << endl << "*****" << endl;
		showHelp = true;		
	}

	if (!showHelp) {

		BedIntersectPE *bi = new BedIntersectPE(bedAFile, bedBFile, overlapFraction, 
												searchType, forceStrand, inputIsBam, outputIsBam);
		bi->DetermineBedPEInput();
		return 0;
	}
	else {
		ShowHelp();
	}
}


void ShowHelp(void) {

	cerr << endl << "Program: " << PROGRAM_NAME << " (v" << VERSION << ")" << endl;
	
	cerr << "Author:  Aaron Quinlan (aaronquinlan@gmail.com)" << endl;

	cerr << "Summary: Report overlaps between a BEDPE file and a BED file." << endl << endl;

	cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -a <BEDPE> -b <BED>" << endl << endl;

	cerr << "Options: " << endl;

	cerr << "\t-abam\t"			<< "The A input file is in BAM format.  Output will be BAM as well." << endl << endl;

	cerr << "\t-outbam\t"		<< "Write output as BAM. This will report those BAM alignments." << endl;
	cerr 						<< "\t\tthat meet the intersect criteria.  BAM entry is written once." << endl << endl;
	
	cerr << "\t-f\t"	    			<< "Minimum overlap required as fraction of A (e.g. 0.05)." << endl;
	cerr 								<< "\t\tDefault is 1E-9 (effectively 1bp)." << endl << endl;

	cerr << "\t-s\t"	    			<< "Enforce strandedness when finding overlaps." << endl;
	cerr 								<< "\t\tDefault is to ignore stand." << endl;
	cerr 								<< "\t\tNot applicable with -type inspan or -type outspan." << endl << endl;
	
	cerr << "\t-type \t"				<< "Approach to reporting overlaps between BEDPE and BED." << endl << endl;
	cerr 								<< "\t\teither\tReport overlaps if either end of A overlaps B." << endl;
	cerr									<< "\t\t\t- Default." << endl;
	cerr 								<< "\t\tneither\tReport A if neither end of A overlaps B." << endl;
	cerr 	 							<< "\t\tboth\tReport overlaps if both ends of A overlap  B." << endl;
	cerr								<< "\t\txor\tReport overlaps if one and only one end of A overlaps B." << endl;
	cerr								<< "\t\tispan\tReport overlaps between [end1, start2] of A and B." << endl; 
	cerr									<< "\t\t\t- Note: If chrom1 <> chrom2, entry is ignored." << endl;
	cerr								<< "\t\tospan\tReport overlaps between [start1, end2] of A and B." << endl; 
	cerr									<< "\t\t\t- Note: If chrom1 <> chrom2, entry is ignored." << endl;
	cerr								<< "\t\tnotispan\tReport A if ispan of A doesn't overlap B." << endl; 
	cerr									<< "\t\t\t\t- Note: If chrom1 <> chrom2, entry is ignored." << endl;
	cerr								<< "\t\tnotospan\tReport A if ospan of A doesn't overlap B." << endl; 
	cerr									<< "\t\t\t\t- Note: If chrom1 <> chrom2, entry is ignored." << endl << endl;
		
	cerr << "Refer to the BEDTools manual for BEDPE format." << endl << endl;

	exit(1);
}
