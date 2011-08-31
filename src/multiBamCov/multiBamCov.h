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
                bool keepDuplicates, bool keepFailedQC);

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
    

    map<string, int> bamFileMap;
    
    void LoadBamFileMap(void);
    void ReportCounts(const vector<int> &counts);
};

#endif /* MULTIBAMCOV_H */
