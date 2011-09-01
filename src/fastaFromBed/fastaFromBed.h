/*****************************************************************************
  fastaFromBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef FASTAFROMBED_H
#define FASTAFROMBED_H

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
class Bed2Fa {

public:

    // constructor
    Bed2Fa(bool useName, const string &dbFile, const string &bedFile, const string &fastaOutFile,
        bool useFasta, bool useStrand);

    // destructor
    ~Bed2Fa(void);

    void ExtractDNA();
    void ReportDNA(const BED &bed, string &dna);


private:

    bool _useName;
    string _dbFile;
    string _bedFile;
    string _fastaOutFile;
    bool _useFasta;
    bool _useStrand;

    // instance of a bed file class.
    BedFile  *_bed;
    ostream *_faOut;
};

#endif
