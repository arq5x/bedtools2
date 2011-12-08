/*****************************************************************************
  sortBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;


//************************************************
// Class methods and elements
//************************************************
class BedSort {

public:

    // constructor
    BedSort(string &bedFile, bool printHeader);

    // destructor
    ~BedSort(void);

    void SortBed();             // the default.  sorts by chrom (asc.) then by start (asc.)
    void SortBedBySizeAsc();
    void SortBedBySizeDesc();
    void SortBedByChromThenSizeAsc();
    void SortBedByChromThenSizeDesc();
    void SortBedByChromThenScoreAsc();
    void SortBedByChromThenScoreDesc();

private:
    string _bedFile;

    // instance of a bed file class.
    BedFile *_bed;
};
