/*****************************************************************************
  multiBamCov.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef MULTICOVBAM_H
#define MULTICOVBAM_H

#include "bedFile.h"
#include "api/BamMultiReader.h"
#include "BlockedIntervals.h"
using namespace BamTools;


#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class MultiCovBam {

public:

    // constructor
    MultiCovBam(const vector<string> &bam_files, const string bed_file, 
                int minQual, bool properOnly, 
                bool keepDuplicates, bool keepFailedQC,
                bool obeySplits, bool sameStrand,
                bool diffStrand, float overlapFraction,
                bool reciprocal);

    // destructor
    ~MultiCovBam(void);

    void CollectCoverage();

private:

    //------------------------------------------------
    // private attributes
    //------------------------------------------------
    vector<string> _bam_files;
    string _bed_file;
	BedFile *_bed;
	
	// attributes to control what is counted
    int _minQual;
    bool _properOnly;
    bool _keepDuplicates;
    bool _keepFailedQC;
    bool _obeySplits;
    bool _sameStrand;
    bool _diffStrand;
    float _overlapFraction;
    bool _reciprocal;


    map<string, int> bamFileMap;
    bool FindBlockedOverlaps(const BED &a, const vector<BED> &a_blocks, 
        const BED &hit);
    void LoadBamFileMap(void);
    void ReportCounts(const vector<int> &counts);
};

#endif /* MULTIBAMCOV_H */
