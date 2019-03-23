/*****************************************************************************
  maskFastaFromBed.h
  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef MASKFASTAFROMBED_H
#define MASKFASTAFROMBED_H

#include "bedFile.h"
#include "sequenceUtils.h"
#include "split.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cctype>   /* for tolower */

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class MaskFastaFromBed {

public:

    // constructor
    MaskFastaFromBed(const string &fastaInFile,  const string &bedFile,
                     const string &fastaOutFile, bool softMask, char maskChar,
                     bool useFullHeader);

    // destructor
    ~MaskFastaFromBed(void);


private:

    bool _softMask;
    bool _useFullHeader;

    string _fastaInFile;
    string _bedFile;
    string _fastaOutFile;
    char   _maskChar;     // typically "N", but user's can choose something else, e.g., "X"

    // instance of a bed file class.
    BedFile *_bed;

    void MaskFasta();

    void PrettyPrintChrom(ofstream &out, string chrom, const string &sequence, CHRPOS width);

};

#endif /* MASKFASTAFROMBED */
