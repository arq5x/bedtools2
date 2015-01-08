/*****************************************************************************
  poolBed.h

  (c) 2015 - Pierre Lindenbaum PhD
  @yokofakun http://plindenbaum.blogspot.com
  Univ. Nantes, France

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <algorithm>  // for binary search
using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedPool {

public:

    // constructor
    BedPool(string &bedFile, string& outfileprefix,size_t num_split);

    // destructor
    ~BedPool(void);

private:

    string _bedFile;
    string _outfileprefix;
    size_t _num_split;

    // The BED file from which to compute coverage.
    BedFile *_bed;
    void doWork();
};

