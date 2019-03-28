/*****************************************************************************
  reldist.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef RELDIST_H
#define RELDIST_H

#include "bedFile.h"
#include "chromsweep.h"
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "BlockedIntervals.h"
#include "BamAncillary.h"
using namespace BamTools;


#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
using namespace std;



class RelativeDistance {

public:

    // constructor
    RelativeDistance(string bedAFile, 
                     string bedBFile,
                     bool _summary);

    // destructor
    ~RelativeDistance(void);

private:

    //------------------------------------------------
    // private attributes
    //------------------------------------------------
    string _bedAFile;
    string _bedBFile;
    
    map<string, vector<CHRPOS> > _db_midpoints;
    map<double, size_t> _reldists;
    size_t _tot_queries;
    
    // instance of a bed file class.
    BedFile *_bedA, *_bedB;
    bool _summary;

    //------------------------------------------------
    // private methods
    //------------------------------------------------
    void LoadMidpoints();
    void CalculateRelativeDistance();
    void UpdateDistanceSummary(double rel_dist);
    void ReportDistanceSummary();



};

#endif /* RELDIST_H */
