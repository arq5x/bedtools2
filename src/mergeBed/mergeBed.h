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
class BedMerge {

public:

  // constructor 
  BedMerge(string &, bool &, int&);

  // destructor
  ~BedMerge(void);

  void MergeBed();

private:
	
	string bedFile;
	bool numEntries;
	int maxDistance;
	// instance of a bed file class.
	BedFile *bed;

};
