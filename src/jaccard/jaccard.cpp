/*****************************************************************************
  jaccard.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "jaccard.h"

/************************************
Helper functions
************************************/
size_t Jaccard::GetTotalIntersection(const BED &a, const vector<BED> &hits) 
{

    size_t intersection = 0;
    
    vector<BED>::const_iterator h       = hits.begin();
    vector<BED>::const_iterator hitsEnd = hits.end();
    for (; h != hitsEnd; ++h) {
        CHRPOS s = max(a.start, h->start);
        CHRPOS e = min(a.end, h->end);
        intersection += (e - s);
    }
    return intersection;
}

/*
    Constructor
*/
Jaccard::Jaccard(string bedAFile, string bedBFile, 
                 float overlapFraction, bool reciprocal)
{
    _bedAFile            = bedAFile;
    _bedBFile            = bedBFile;
    _overlapFraction     = overlapFraction;
    _reciprocal          = reciprocal;

        
    CalculateJaccard();
}


/*
    Destructor
*/
Jaccard::~Jaccard(void) {
}


unsigned long Jaccard::GetIntersection(size_t &n_intersections) {
    
    _bedA = new BedFile(_bedAFile);
    _bedB = new BedFile(_bedBFile);
    
    unsigned long I = 0;

    ChromSweep sweep = ChromSweep(_bedA, _bedB, 
                                  false, false, 
                                  _overlapFraction, _reciprocal,
                                  true, false);


    pair<BED, vector<BED> > hit_set;
    hit_set.second.reserve(10000);
    while (sweep.Next(hit_set)) {
        I += GetTotalIntersection(hit_set.first, hit_set.second);
        n_intersections += hit_set.second.size();
    }
    return I;
}

void Jaccard::CalculateJaccard() {

    size_t n_intersections = 0;
    unsigned long I = GetIntersection(n_intersections);
    
    unsigned long U = _bedA->getTotalFlattenedLength() + \
                      _bedB->getTotalFlattenedLength();
    
    // header
    cout << "intersection\t"
         << "union\t"
         << "jaccard\t"
         << "n_intersections"
         << endl;
    
    // result
    cout << I << "\t" 
         << U - I << "\t"
         << (float) I / ((float) U - (float) I) << "\t"
         << n_intersections
         << endl;
}


