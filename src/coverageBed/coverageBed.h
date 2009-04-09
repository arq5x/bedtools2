#ifndef GRAPHBED_H
#define GRAPHBED_H

#include "bedFile.h"
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <fstream>

using namespace std;


//***********************************************
// Typedefs
//***********************************************
typedef list<BED> bedList;


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
#endif /* GRAPHBED_H */
