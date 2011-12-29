/*****************************************************************************
  windowMaker.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "genomeFile.h"

using namespace std;


//************************************************
// Class methods and elements
//************************************************
class WindowMaker {

public:

  // constructor
  WindowMaker(string &genomeFile, uint32_t size, uint32_t step);

  // destructor
  ~WindowMaker(void);

  void MakeWindows();

private:
    string _genomeFile;
    GenomeFile *_genome;
    uint32_t _size; 
    uint32_t _step;
};
