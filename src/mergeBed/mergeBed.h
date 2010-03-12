/*****************************************************************************
  mergeBed.h

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

void ReportMergedNames(const map<string, bool> &names);

//************************************************
// Class methods and elements
//************************************************
class BedMerge {

public:

  // constructor 
  BedMerge(string &, bool &, int&, bool &, bool &);

  // destructor
  ~BedMerge(void);

  void MergeBed();
  void MergeBedStranded();

private:
	
	string bedFile;
	bool numEntries;
	bool forceStrand;
	bool reportNames;
	int maxDistance;
	// instance of a bed file class.
	BedFile *bed;

};
