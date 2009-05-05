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
