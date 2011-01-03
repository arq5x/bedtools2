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
    MaskFastaFromBed(string &fastaInFile, string &bedFile, string &fastaOutFile, bool &softMask);

    // destructor
    ~MaskFastaFromBed(void);


private:

    bool _softMask;

    string _fastaInFile;
    string _bedFile;
    string _fastaOutFile;

    // instance of a bed file class.
    BedFile *_bed;

    void MaskFasta();

    void PrettyPrintChrom(ofstream &out, string chrom, const string &sequence, int width);

};

#endif /* MASKFASTAFROMBED */
