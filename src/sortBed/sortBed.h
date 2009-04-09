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
class BedSort {

public:

  // constructor 
  BedSort(string &);

  // destructor
  ~BedSort(void);

  void SortBed();

private:
	
	string bedFile;

	// instance of a bed file class.
	BedFile *bed;

};
