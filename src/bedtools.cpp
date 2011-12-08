/*****************************************************************************
  bedtools.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include <vector>
#include <algorithm>    // for swap()
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;


// define our program name
#define PROGRAM_NAME "bedtools"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


int intersect_main(int argc, char* argv[]);
int coverage_main(int argc, char* argv[]);
int merge_main(int argc, char* argv[]);
int substract_main(int argc, char* argv[]);

static int help()
{
	cerr << "\n";
	cerr << "Program: " << PROGRAM_NAME << " (Tools for genome arithmetic.)\n";
	cerr << "Version: " << "2.15.5" << "\n";
	cerr << "Author:  " << "Aaron Quinlan (aaronquinlan@gmail.com)" << "\n\n";
	cerr << "Usage:   bedtools <tool> [options]\n";

	cerr << "\n-Genome arithmetic tools:\n";	
	cerr << "   intersect     Find overlapping/intersecting intervals b/w two files.\n";
	cerr << "   coverage      Compute the coverage of one set of intervals over another.\n";
	cerr << "   merge         Combine overlapping or nearby intervals into a single interval.\n";
	cerr << "   subtract      Remove intervals based on overlaps b/w two files.\n";
	
	cerr << "\n-Format conversion tools:\n";	
	cerr << "   intersect     Find overlapping/intersecting intervals b/w two files.\n";
	
	cerr << "\n-Fasta tools:\n";	
	cerr << "   intersect     Find overlapping/intersecting intervals b/w two files.\n";
	
	cerr << "\n";
	return 1;
}

int main(int argc, char *argv[])
{
    // make sure the user at least entered a sub_command
	if (argc < 2) return help();

    // spawn the proper tool.
	if (strcmp(argv[1], "intersect") == 0)     return intersect_main(argc-1, argv+1);
	else if (strcmp(argv[1], "coverage") == 0) return coverage_main(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0)    return merge_main(argc-1, argv+1);
	else if (strcmp(argv[1], "subtract") == 0) return subtract_main(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;	
}
