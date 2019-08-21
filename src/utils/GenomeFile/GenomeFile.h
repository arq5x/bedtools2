/*****************************************************************************
  GenomeFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef GENOMEFILE_H
#define GENOMEFILE_H

#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <algorithm> // for bsearch lower_bound()
#include "api/BamReader.h"
#include "api/BamAux.h"
using namespace BamTools;

using namespace std;

typedef int64_t CHRPOS;

// typedef for mapping b/w chrom name and it's size in b.p.
typedef map<string, CHRPOS, std::less<string> > chromToSizes;


class GenomeFile {

public:

    // Constructor using a file
    GenomeFile(const string &genomeFile);
    
    // Constructor using a vector of BamTools RefVector
    GenomeFile(const RefVector &genome);

    // Destructor
    ~GenomeFile(void);

    // load a GENOME file into a map keyed by chrom. value is size of chrom.
    void loadGenomeFileIntoMap();
    
    pair<string, CHRPOS> projectOnGenome(CHRPOS genome_pos);
    
    uint64_t getChromSize(const string &chrom);     // return the size of a chromosome
    uint64_t getGenomeSize(void);              // return the total size of the geonome
    vector<string> getChromList();             // return a list of chrom names
    int getNumberOfChroms();                   // return the number of chroms
    string getGenomeFileName();                // return the name of the genome file




private:
    string  _genomeFile;
    chromToSizes _chromSizes;
    vector<string> _chromList;

    // projecting chroms into a single coordinate system
    CHRPOS _genomeLength;
    vector<CHRPOS> _startOffsets;
    
};

#endif /* GENOMEFILE_H */
