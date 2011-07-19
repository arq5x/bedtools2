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
    Bed2Fa(const string &dbFile, const string &bedFile, 
           bool useFasta, bool useStrand, bool useName, bool useFull);

    // destructor
    ~Bed2Fa(void);

    void ExtractDNA();
    void ReportDNA(const BED &bed, string &dna);


private:
    string _dbFile;
    string _bedFile;
    bool _useFasta;
    bool _useStrand;  // should we extract a specific strand?
    bool _useName;    // should we use the BED name for the FASTA header?
    bool _useFull;    // should we use the full BED entry for the tabular output?

    // instance of a bed file class.
    BedFile  *_bed;
};

#endif
