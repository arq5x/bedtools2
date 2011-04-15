/*****************************************************************************
  nucBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef NUCBED_H
#define NUCBED_H

#include "bedFile.h"
#include "sequenceUtils.h"
#include "Fasta.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class NucBed {

public:

    // constructor
    NucBed(string &dbFile, string &bedFile, bool printSeq);

    // destructor
    ~NucBed(void);

    void ProfileDNA();


private:
    string _dbFile;
    string _bedFile;
    bool   _printSeq;

    // instance of a bed file class.
    BedFile  *_bed;
    
    void ReportDnaProfile(const BED& bed, const string &sequence, int seqLength);
};

#endif
