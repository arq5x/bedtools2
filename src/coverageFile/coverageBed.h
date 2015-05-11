/*****************************************************************************
  coverageBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef COVERAGEBED_H
#define COVERAGEBED_H

#include "bedFile.h"

#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BlockedIntervals.h"
using namespace BamTools;

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

//************************************************
// Class methods and elements
//************************************************
class BedCoverage {

public:

    // constructor
    BedCoverage(string &bedAFile, string &bedBFile, bool sameStrand, bool diffStrand, bool writeHistogram,
                bool bamInput, bool obeySplits, bool eachBase, bool countsOnly);

    // destructor
    ~BedCoverage(void);

private:

    // input files.
    string _bedAFile;
    string _bedBFile;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    // do we care about same or opposite strandedness when counting coverage?
    bool _sameStrand;
    bool _diffStrand;

    // should we write a histogram for each feature in B?
    bool _writeHistogram;

    // are we dealing with BAM input for "A"?
    bool _bamInput;

    // should we split BED/BAM into discrete blocks?
    bool _obeySplits;

    // should discrete coverage be reported for each base in each feature?
    bool _eachBase;
    
    // should we just count overlaps and not try to describe the breadth?
    bool _countsOnly;

    // private function for reporting coverage information
    void ReportCoverage();
    
    // private function for reporting overlap counts
    void ReportCounts();

    void CollectCoverageBed();

    void CollectCoverageBam(string bamFile);
};
#endif /* COVERAGEBED_H */
