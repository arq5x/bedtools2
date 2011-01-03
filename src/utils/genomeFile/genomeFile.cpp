/*****************************************************************************
  genomeFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "genomeFile.h"


GenomeFile::GenomeFile(const string &genomeFile) {
    _genomeFile = genomeFile;
    loadGenomeFileIntoMap();
}

// Destructor
GenomeFile::~GenomeFile(void) {
}


void GenomeFile::loadGenomeFileIntoMap() {

    string genomeLine;
    int lineNum = 0;
    vector<string> genomeFields;            // vector for a GENOME entry

    // open the GENOME file for reading
    ifstream genome(_genomeFile.c_str(), ios::in);
    if ( !genome ) {
        cerr << "Error: The requested genome file (" << _genomeFile << ") could not be opened. Exiting!" << endl;
        exit (1);
    }

    while (getline(genome, genomeLine)) {

        Tokenize(genomeLine,genomeFields);  // load the fields into the vector
        lineNum++;

        // ignore a blank line
        if (genomeFields.size() > 0) {
            if (genomeFields[0].find("#") == string::npos) {

                // we need at least 2 columns
                if (genomeFields.size() >= 2) {
                    char *p2End;
                    long c2;
                    // make sure the second column is numeric.
                    c2 = strtol(genomeFields[1].c_str(), &p2End, 10);

                    // strtol  will set p2End to the start of the string if non-integral, base 10
                    if (p2End != genomeFields[1].c_str()) {
                        string chrom       = genomeFields[0];
                        int size           = atoi(genomeFields[1].c_str());
                        _chromSizes[chrom] = size;
                        _chromList.push_back(chrom);
                    }
                }
                else {
                    cerr << "Less than the req'd two fields were encountered in the genome file (" << _genomeFile << ")";
                    cerr << " at line " << lineNum << ".  Exiting." << endl;
                    exit (1);
                }
            }
        }
        genomeFields.clear();
    }
}


int GenomeFile::getChromSize(const string &chrom) {
    chromToSizes::const_iterator chromIt = _chromSizes.find(chrom);
    if (chromIt != _chromSizes.end())
        return _chromSizes[chrom];
    else
        return -1;  // chrom not found.
}

vector<string> GenomeFile::getChromList() {
    return _chromList;
}

int GenomeFile::getNumberOfChroms() {
    return _chromList.size();
}

string GenomeFile::getGenomeFileName() {
    return _genomeFile;
}
