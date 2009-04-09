#include "bedFile.h"
#include <vector>
#include <algorithm>
#include <map>
#include <list>
#include <iostream>
#include <fstream>
#include <limits>

using namespace std;


//***********************************************
// Typedefs
//***********************************************
typedef list<BED> bedList;


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
