#include "slopBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "slopBed"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void ShowHelp(void);

int main(int argc, char* argv[]) {

	// our configuration variables
	bool showHelp = false;

	// input files
	string bedFile;
	string genomeFile;
	
	bool haveBed = false;
	bool haveGenome = false;
	bool haveLeft = false;
	bool haveRight = false;
	bool haveBoth = false;	
	
	bool forceStrand = false;
	int leftSlop = 0;
	int rightSlop = 0;

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
		else if(PARAMETER_CHECK("-g", 2, parameterLength)) {
			haveGenome = true;
			genomeFile = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-l", 2, parameterLength)) {
			haveLeft = true;
			leftSlop = atoi(argv[i + 1]);
			i++;
		}
		else if(PARAMETER_CHECK("-r", 2, parameterLength)) {
			haveRight = true;			
			rightSlop = atoi(argv[i + 1]);
			i++;
		}
		else if(PARAMETER_CHECK("-b", 2, parameterLength)) {
			haveBoth = true;
			leftSlop = atoi(argv[i + 1]);
			rightSlop = atoi(argv[i + 1]);			
			i++;
		}				
		else if(PARAMETER_CHECK("-s", 2, parameterLength)) {
			forceStrand = true;
			i++;
		}	
		else {
		  cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
			showHelp = true;
		}		
	}

	// make sure we have both input files
	if (!haveBed || !haveGenome) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need both a BED (-i) and a genome (-g) file. " << endl << "*****" << endl;
	  showHelp = true;
	}
	if (!haveLeft && !haveRight && !haveBoth) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need -l and -r together or -b alone. " << endl << "*****" << endl;
	  showHelp = true;
	}
	if ((!haveLeft && haveRight) || (haveLeft && !haveRight)) {
	  cerr << endl << "*****" << endl << "*****ERROR: Need both -l and -r. " << endl << "*****" << endl;
	  showHelp = true;
	}	
	if (forceStrand && (!(haveLeft) || !(haveRight))) {
	  cerr << endl << "*****" << endl << "*****ERROR: Must supply -l and -r with -s. " << endl << "*****" << endl;
	  showHelp = true;	
	}	
	if (!showHelp) {
		BedSlop *bc = new BedSlop(bedFile, genomeFile, forceStrand, leftSlop, rightSlop);
		bc->ProcessBed();		
		
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
	cerr << "Description: Add requested base pairs of \"slop\" to each BED entry." << endl << endl;

	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -g <genome> -i <bed>" << endl << endl;
	
	cerr << "OPTIONS: " << endl;
	cerr << "\t" << "-l\t\t\t"            	<< "The number of base pairs to subtract from the start coordinate." << endl;
	cerr << "\t" << "-r\t\t\t"            	<< "The number of base pairs to add to the end coordinate." << endl;
	cerr << "\t" << "-b\t\t\t"            	<< "Increase the BED entry by the same number base pairs in each direction." << endl;
	cerr << "\t" << "-s\t\t\t"            	<< "Define -l and -r based on strand.  E.g. if used, -l 500 for a negative-stranded " << endl << "\t\t\t\tfeature, it will add 500 bp downstream.  Default = false." << endl << endl;	

	cerr << "NOTES: " << endl;
	cerr << "\t1.  Starts will be corrected to 0 if the requested slop would force it below 0." << endl;
	cerr << "\t2.  Ends will be corrected to the chromosome length if the requested slop would force it above the max chrom length." << endl << endl;
	cerr << "\t3.  The genome file should tab delimited and structured as follows: <chr><TAB><size>. For example, Mus musculus:" << endl;
	cerr << "\t\tchr1\t197195432" << endl;
	cerr << "\t\tchr2\t181748087" << endl;
	cerr << "\t\t..." << endl;
	cerr << "\t\tchrY_random\t58682461" << endl << endl;
	cerr << "\t4.  ***Only tab-delimited BED3 - BED6 formats allowed.***"<< endl << endl;

	cerr << "TIPS:" << endl;
	cerr << "\tOne can use the UCSC Genome Browser's MySQL database to extract chromosome sizes. For example, H. sapiens:" << endl << endl;
	cerr << "\tmysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \"select chrom, size from hg18.chromInfo\"  > hg18.genome" 
		<< endl << endl;
		
	
	// end the program here
	exit(1);
}
