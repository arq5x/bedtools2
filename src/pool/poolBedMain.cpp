/*****************************************************************************
  dealBed.h

  (c) 2015 - Pierre Lindenbaum PhD
  @yokofakun http://plindenbaum.blogspot.com
  Univ. Nantes, France

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "poolBed.h"
#include "version.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "bedtools pool"


// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
void pool_help(void);

int pool_main(int argc, char* argv[]) {

    // our configuration variables
    bool showHelp = false;
    bool haveBed  = false;
    // input files
    string bedFile("stdin");
    string outfileprefix("_pool");
    size_t num_split(10UL);


    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);
	if(PARAMETER_CHECK("-i", 2, parameterLength)) {
		    if ((i+1) < argc) {
		        bedFile = argv[i + 1];
		        haveBed=true;
		        i++;
		    }
		}
	else if(PARAMETER_CHECK("-p", 2, parameterLength)) {
		    if ((i+1) < argc) {
		        outfileprefix.assign(argv[i + 1]);
		        i++;
		    }
		}
	else if(PARAMETER_CHECK("-n", 2, parameterLength)) {
		    if ((i+1) < argc) {
		        num_split = atof(argv[i + 1]);
		        i++;
		    }
		}
        else if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }
    
    if(outfileprefix.empty())
        {
        fprintf(stderr,"output prefix is empty\n");
        exit(EXIT_FAILURE);
        }
    
    if(showHelp) pool_help();
    
    if( num_split < 0UL) num_split = 1UL;
        
    if (!showHelp) {
        BedPool *bc = new BedPool(bedFile,
        	outfileprefix,num_split);
        
        delete bc;
        return 0;
    }
    else {
        pool_help();
    }
    return 0;
}

void pool_help(void) {
    cerr << "\nTool:    bedtools pool" << endl;
    cerr << "Version: " << VERSION << "\n";    
    cerr << "Summary: Groups items of a Bed file in " << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bed> -n number-of-files" << endl << endl;

    cerr << "Options: " << endl;
    cerr << "\t-n\t"                << "Number of pools to create." << endl;
    cerr << "\t-p\t"                << "output file prefix" << endl;

    // end the program here
    exit(1);
}
