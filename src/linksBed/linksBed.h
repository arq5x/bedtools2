/*****************************************************************************
  linksBed.h

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
class BedLinks {

public:

    // constructor
    BedLinks(string &bedFile, string &base, string &org, string &db);

    // destructor
    ~BedLinks(void);

private:
    string _bedFile;
    string _base;
    string _org;
    string _db;

    // instance of a bed file class.
    BedFile *_bed;

    void WriteURL(BED &bed, string &base);
    void CreateLinks();             // the default.  sorts by chrom (asc.) then by start (asc.)
};
