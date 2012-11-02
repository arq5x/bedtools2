/*****************************************************************************
  chromsweepBed.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef CHROMSWEEP_H
#define CHROMSWEEP_H

#include "bedFile.h"
#include <vector>
#include <list>
#include <queue>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class ChromSweep {

// public interface.
public:

    // A is the query and B is the database
    
    // constructor using existing BedFile pointers
    ChromSweep(BedFile *query, BedFile *db, 
               bool sameStrand = false, bool diffStrand = false, 
               float overlapFraction = 0.0, bool reciprocal = false,
               bool useMergedIntervals = false, bool printHeader = false);
    
    // constructor using filenames
    ChromSweep(string &queryFile, string &dbFile);
    
    // destructor
    ~ChromSweep(void);
    
    // loads next (a pair) with the current query and it's overlaps
    //   next.first is the current query interval
    //   next.second is a vector of the current query's hits.
    // returns true if overlap
    bool Next(pair<BED, vector<BED> > &next);
    
    // Usage:
    //     ChromSweep sweep = ChromSweep(_bedA, _bedB);
    //     pair<BED, vector<BED> > hit_set;
    //     while (sweep.Next(hit_set)) 
    //     {
    //        // magic happens here!
    //        processHits(hit_set.first, hit_set.second);
    //     }
    
// private variables.
private:

    // instances of a bed file class.
    BedFile *_query, *_db;
    float _overlapFraction;
    // do we care about strandedness.
    bool _sameStrand;
    bool _diffStrand;
    // do we care about reciprocal overlap?
    bool _reciprocal;
    // should we merge overlapping intervals before computing overlaps?
    bool _useMergedIntervals;

    /* 
       a cache of still active features from the database file\
       2012-Oct-29: prefer LIST over VECTOR as, when _cache is large,
       the overhead of deleting from a vector is huge.  Deleting from
       the front of a LIST is cheap. Thanks to Neil Kindlon (Quinlan lab)
       for the very important fix. 
    */
    list<BED> _cache;
    // the set of hits in the database for the current query
    vector<BED> _hits;
    // a queue from which we retrieve overlap results.  used by Next()
    queue< pair<BED, vector<BED> > > _results;
    BED _nullBed;
    // an empty BED vector for returning no hits for a given query
    vector<BED> _no_hits;
    // the current query and db features.
    BED _curr_qy, _curr_db;
    // a cache of the current chrom from the query. used to handle chrom changes.
    string _curr_chrom;

// private methods.
private:
    
    void ScanCache();   
    bool ChromChange();
    bool NextQuery();
    bool NextDatabase();
    bool IsValidHit(const BED &query, const BED &db);
};

#endif /* CHROMSWEEP_H */
