#include "bedFile.h"
#include <vector>
#include <algorithm>
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
class BedLinks {

public:

	// constructor 
	BedLinks(string &, string &, string &, string &);

	// destructor
	~BedLinks(void);

	void WriteURL(BED &, string &);
	void LinksBed();				// the default.  sorts by chrom (asc.) then by start (asc.)

private:	
	string bedFile;
	string base;
	string org;
	string db;

	// instance of a bed file class.
	BedFile *bed;

};
