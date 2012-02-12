/*****************************************************************************
  mergeBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits.h>
#include <stdlib.h>
#include "VectorOps.h"

using namespace std;

const int PRECISION = 21;

//************************************************
// Class methods and elements
//************************************************
class BedMerge {

public:

  // constructor
  BedMerge(string &bedFile, bool numEntries, 
           int maxDistance, bool forceStrand, 
           bool reportNames, bool reportScores, const string &scoreOp);

  // destructor
  ~BedMerge(void);

  void MergeBed();
  void MergeBedStranded();

private:

    string _bedFile;
    bool   _numEntries;
    bool   _forceStrand;
    bool   _reportNames;
    bool   _reportScores;
    string _scoreOp;
    int    _maxDistance;
    // instance of a bed file class.
    BedFile *_bed;

    void Report(string chrom, int start, int end, const vector<string> &names, const vector<string> &scores, int mergeCount);
    void ReportStranded(string chrom, int start, int end, const vector<string> &names, const vector<string> &scores, int mergeCount, string strand);
    void ReportMergedNames(const vector<string> &names);
    void ReportMergedScores(const vector<string> &scores);
    
};
