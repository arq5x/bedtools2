#ifndef FASTAFROMBED_H
#define FASTAFROMBED_H

#include "bedFile.h"
#include "sequenceUtils.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class Bed2Fa {

public:
	
	// constructor 
	Bed2Fa(bool &, string &, string &, string &);

	// destructor
	~Bed2Fa(void);

	void ExtractDNA();
	
private:
	
	bool useName;
	string dbFile;
	string bedFile;
	string fastaOutFile;
	
	// instance of a bed file class.
	BedFile *bed;

};

#endif
