/*****************************************************************************
  genomeFile.h

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

using namespace std;


// typedef for mapping b/w chrom name and it's size in b.p.
typedef map<string, int, std::less<string> > chromToSizes;


class GenomeFile {

public:

    // Constructor
    GenomeFile(const string &genomeFile);

    // Destructor
    ~GenomeFile(void);

    // load a GENOME file into a map keyed by chrom. value is size of chrom.
    void loadGenomeFileIntoMap();

    int getChromSize(const string &chrom);  // return the size of a chromosome
    vector<string> getChromList();          // return a list of chrom names
    int getNumberOfChroms();                // return the number of chroms
    string getGenomeFileName();             // return the name of the genome file



private:
    string  _genomeFile;
    chromToSizes _chromSizes;
    vector<string> _chromList;
};

#endif /* GENOMEFILE_H */
