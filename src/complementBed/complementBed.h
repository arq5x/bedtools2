/*****************************************************************************
  complementBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "bedFile.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits.h>
#include <stdlib.h>

using namespace std;


//************************************************
// Class methods and elements
//************************************************
class BedComplement {

public:

  // constructor 
  BedComplement(string &, string &);

  // destructor
  ~BedComplement(void);

  void ComplementBed();

private:
	
	string bedFile;
	string genomeFile;
		
	// instance of a bed file class.
	BedFile *bed;

};
