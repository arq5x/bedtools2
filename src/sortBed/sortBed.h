#include "bedFile.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;


//************************************************
// Class methods and elements
//************************************************
class BedSort {

public:

	// constructor 
	BedSort(string &);

	// destructor
	~BedSort(void);

	// write BED to stdout
	void reportBed(const BED &);

	void SortBed();				// the default.  sorts by chrom (asc.) then by start (asc.)
	void SortBedBySizeAsc();
	void SortBedBySizeDesc();
	void SortBedByChromThenSizeAsc();
	void SortBedByChromThenSizeDesc();
	void SortBedByChromThenScoreAsc();
	void SortBedByChromThenScoreDesc();
	
private:	
	string bedFile;

	// instance of a bed file class.
	BedFile *bed;

};
