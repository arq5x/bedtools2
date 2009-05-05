#include <iostream>	
#include "linksBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "linksBed"


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

	/* Defaults for Quinlan
	//string org = "mouse";
	//string db = "mm9";
	//string base = "http://mirror.bioch.virginia.edu";
	*/
	
	/* Defaults for everyone else */
	string org = "human";
	string db = "hg18";
	string base = "http://genome.ucsc.edu";
	
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

		if(PARAMETER_CHECK("-i", 2, parameterLength)) {
			haveBed = true;
			bedFile = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-base", 5, parameterLength)) {
			base = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-org", 4, parameterLength)) {
			org = argv[i + 1];
			i++;
		}
		else if(PARAMETER_CHECK("-db", 3, parameterLength)) {
			db = argv[i + 1];
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

	if (!showHelp) {
		BedLinks *bl = new BedLinks(bedFile, base, org, db);
		bl->LinksBed();			
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
	cerr << "Description: Creates HTML links from a BED file." << endl << endl;
	cerr << "Usage: " << PROGRAM_NAME << " [OPTIONS] -i <input.bed> > out.html" << endl << endl;

	cerr << "OPTIONS: " << endl;
	cerr << "\t" << "-base\t\t"	<< "The browser basename.  Default: http://genome.ucsc.edu " << endl << endl;
	cerr << "\t" << "-org\t\t"	<< "The organism. Default: human" << endl << endl;
	cerr << "\t" << "-db\t\t"	<< "The build.  Default: hg18" << endl << endl;

	cerr << "NOTES: " << endl;
	cerr << "\t" << "-i stdin\t"	<< "Allows BED file A to be read from stdin.  E.g.: cat a.bed | linksBed -i stdin > out.html" << endl << endl;
	cerr << "\t***Only BED3 - BED6 formats allowed.***"<< endl << endl;
	
	cerr << "EXAMPLE: " << endl;
	cerr << "\t" << "By default, the links created will point to the human (hg18) UCSC browser.  If you have a local mirror, you can override this behavior";
	cerr << "by supplying the -base, -org, and -db options."  << endl;
	cerr << "\t" << "For example, if the main URL of your local mirror for mouse MM9 is called: http://mymirror.myuniversity.edu, then you would use the following parameters:" << endl;
	cerr << "\t\t" << "-base http://mymirror.myuniversity.edu" << endl;
	cerr << "\t\t" << "-org mouse" << endl; 
	cerr << "\t\t" << "-db mm9" << endl;
	// end the program here
	exit(1);

}
