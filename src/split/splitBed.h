/*****************************************************************************
  splitBed.h

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


//************************************************
// Class methods and elements
//************************************************
class BedSplit {

private:
    std::string bedFileName;
    std::string outfileprefix;
    unsigned int num_chuncks;
    int bedType;
    std::vector<BED> items;
    
    void usage(std::ostream& out);
    std::FILE* saveFileChunk(std::string& name,size_t file_index);
    void saveBedItems(void* data,size_t file_index);
    void loadBed();
    int doSimpleSplit();
    int doEuristicSplitOnTotalSize();

public:
    // constructor
    BedSplit();
    // destructor
    ~BedSplit(void);
    //main method
    int main(int argc,char** argv);
};

