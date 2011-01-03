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

#include "BamReader.h"
#include "BamAux.h"
#include "BamAncillary.h"
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
    BedCoverage(string &bedAFile, string &bedBFile, bool forceStrand, bool writeHistogram,
                bool bamInput, bool obeySplits, bool eachBase);

    // destructor
    ~BedCoverage(void);

private:

    // input files.
    string _bedAFile;
    string _bedBFile;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    // do we care about strandedness when counting coverage?
    bool _forceStrand;

    // should we write a histogram for each feature in B?
    bool _writeHistogram;

    // are we dealing with BAM input for "A"?
    bool _bamInput;

    // should we split BED/BAM into discrete blocks?
    bool _obeySplits;

    // should discrete coverage be reported for each base in each feature?
    bool _eachBase;

    // private function for reporting coverage information
    void ReportCoverage();

    void CollectCoverageBed();

    void CollectCoverageBam(string bamFile);
};
#endif /* COVERAGEBED_H */
