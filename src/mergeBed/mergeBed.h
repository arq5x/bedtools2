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
#include <iostream>
#include <fstream>
#include <limits.h>
#include <stdlib.h>

using namespace std;

void ReportMergedNames(const map<string, bool> &names);

//************************************************
// Class methods and elements
//************************************************
class BedMerge {

public:

  // constructor
  BedMerge(string &bedFile, bool &numEntries, int &maxDistance, bool &forceStrand, bool &reportNames);

  // destructor
  ~BedMerge(void);

  void MergeBed();
  void MergeBedStranded();

private:

    string _bedFile;
    bool _numEntries;
    bool _forceStrand;
    bool _reportNames;
    int _maxDistance;
    // instance of a bed file class.
    BedFile *_bed;

    void Report(string chrom, int start, int end, const map<string, bool> &names, int mergeCount);
    void ReportStranded(string chrom, int start, int end, const map<string, bool> &names, int mergeCount, string strand);
};
