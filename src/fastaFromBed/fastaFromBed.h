/*****************************************************************************
  fastaFromBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
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
	Bed2Fa(bool &, string &, string &, string &, bool &);

	// destructor
	~Bed2Fa(void);

	void ExtractDNA();

private:
	
	bool useName;
	string dbFile;
	string bedFile;
	string fastaOutFile;
	bool useFasta;
	
	// instance of a bed file class.
	BedFile *bed;

};

#endif
