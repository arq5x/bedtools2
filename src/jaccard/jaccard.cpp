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


unsigned long Jaccard::GetUnion() {

    // create new BED file objects for A and B
    _bedA = new BedFile(_bedAFile);
    _bedB = new BedFile(_bedBFile);

    unsigned long U = 0;
    BED bed;    
    _bedA->Open();
    while (_bedA->GetNextMergedBed(bed)) {
        U += bed.end - bed.start;
    }
    
    _bedB->Open();
    while (_bedB->GetNextMergedBed(bed)) {
        U += bed.end - bed.start;
    }
    
    return U;
}

unsigned long Jaccard::GetIntersection() {
    
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
    }
    return I;
}

void Jaccard::CalculateJaccard() {

    unsigned long U = GetUnion();
    delete _bedA;
    delete _bedB;
    unsigned long I = GetIntersection();
    
    // header
    cout << "intersection\t"
         << "union\t"
         << "jaccard"
         << endl;
    
    // result
    cout << I << "\t" 
         << U - I << "\t"
         << (float) I / ((float) U - (float) I)
         << endl;
}


