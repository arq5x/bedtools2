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
