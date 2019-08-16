
/*****************************************************************************
  NewGenomeFile.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef NEW_GENOMEFILE_H
#define NEW_GENOMEFILE_H

#include <algorithm> // for bsearch lower_bound()

#include "BedtoolsTypes.h"

#include "api/BamReader.h"
#include "api/BamAux.h"


class NewGenomeFile {

public:

     NewGenomeFile(const string &genomeFileName);
     NewGenomeFile(const BamTools::RefVector &genome);
    ~NewGenomeFile(void);

    // load a GENOME file into a map keyed by chrom. value is a pair<int, int> of id and size.
    void loadGenomeFileIntoMap();
    
    bool projectOnGenome(CHRPOS genome_pos, string &chrom, CHRPOS &start);
    
    CHRPOS getGenomeSize(void) const { return _genomeLength; }                // return the total size of the genome
    CHRPOS getChromSize(const string &chrom);  // return the size of a chromosome
    CHRPOS getChromSize(const string &chrom) const;  // return the size of a chromosome
    uint32_t getChromId(const string &chrom); // return chromosome's sort order
    const vector<string> &getChromList() const { return _chromList; }  // return a list of chrom names
    uint32_t getNumberOfChroms() const { return (uint32_t)(_chromList.size() -1); }//the -1 excludes the blank chrom added for unmapped reads
    const string &getGenomeFileName() const { return _genomeFileName; }
    bool hasChrom(const string &chrom) const { return _chromSizeIds.find(chrom) != _chromSizeIds.end(); }




private:
    string  _genomeFileName;
    istream   *_genomeFile;
    typedef map<string, pair<CHRPOS, int> > lookupType;
    lookupType _chromSizeIds;
    vector<string> _chromList;
    int _maxId;

    // projecting chroms onto a single coordinate system
    CHRPOS _genomeLength;
    vector<CHRPOS> _startOffsets;
    
    //cache members for quick lookup
    string _currChromName;
    CHRPOS _currChromSize;
    int _currChromId;

};

#endif /* GENOMEFILE_H */
