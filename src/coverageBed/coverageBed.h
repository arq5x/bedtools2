#ifndef	COVERAGEBED_H
#define COVERAGEBED_H

#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedGraph {

public:

	// constructor 
	BedGraph(string &, string &);

	// destructor
	~BedGraph(void);
	
	void GraphBed();
	
	
private:

	string bedAFile;
	string bedBFile;

	// instance of a bed file class.
	BedFile *bedA, *bedB;

};
#endif /* COVERAGEBED_H */
