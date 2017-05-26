/*****************************************************************************
  clusterBed.h

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
#include <fstream>
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;


//************************************************
// Class methods and elements
//************************************************
class BedCluster {

public:

  // constructor
  BedCluster(string &bedFile, int maxDistance, bool forceStrand);
  // destructor
  ~BedCluster(void);
  // find clusters
  void ClusterBed();
  // find clusters based on strand
  void ClusterBedStranded();

private:
    string _bedFile;
    bool   _forceStrand;
    int    _maxDistance;
    // instance of a bed file class.
    BedFile *_bed;    
};
