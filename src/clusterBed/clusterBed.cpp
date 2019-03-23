/*****************************************************************************
  clusterBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "clusterBed.h"

// = Constructor =
BedCluster::BedCluster(string &bedFile, 
                   int  maxDistance, 
                   bool forceStrand) 
    :
    _bedFile(bedFile),
    _forceStrand(forceStrand),
    _maxDistance(maxDistance)
{
    _bed = new BedFile(bedFile);
    if (_forceStrand == false)
        ClusterBed();
    else
        ClusterBedStranded();
}


// = Destructor =
BedCluster::~BedCluster(void) 
{}


// = Cluster overlapping or nearby BED entries =
void BedCluster::ClusterBed() {

    uint32_t cluster_id = 0;
    BED prev, curr;
    CHRPOS end   = -1;
    
    _bed->Open();
    while (_bed->GetNextBed(curr, true)) { // true = force sorted intervals
        if (_bed->_status != BED_VALID)
            continue;            

        CHRPOS distance = (curr.start - end);
        
        // new cluster, no overlap
        if ( (distance > _maxDistance) || (curr.chrom != prev.chrom) ) 
        {
            cluster_id++;
            end   = curr.end;
        }
        else {
            if (curr.end > end)
                end   = curr.end;
        }
        prev = curr;
        _bed->reportBedTab(curr);
        printf("%d\n", cluster_id);
    }
}


// = Cluster overlapping BED entries, accounting for strandedness =
void BedCluster::ClusterBedStranded() {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bed->loadBedFileIntoMapNoBin();
    
    uint32_t cluster_id = 0;

    // loop through each chromosome and merge their BED entries
    masterBedMapNoBin::const_iterator m    = _bed->bedMapNoBin.begin();
    masterBedMapNoBin::const_iterator mEnd = _bed->bedMapNoBin.end();
    for (; m != mEnd; ++m) {
        
        // bedList is already sorted by start position.
        string chrom        = m->first;
        vector<BED> bedList = m->second;
    
        // make a list of the two strands to merge separately.
        vector<string> strands(2);
        strands[0] = "+";
        strands[1] = "-";
        // do two passes, one for each strand.
        for (unsigned int s = 0; s < strands.size(); s++) {
            // cluster overlapping features for this chromosome.
            CHRPOS end   = -1;
            BED prev;    
            vector<BED>::const_iterator bedItr = bedList.begin();
            vector<BED>::const_iterator bedEnd = bedList.end();
            for (; bedItr != bedEnd; ++bedItr) {
                // if forcing strandedness, move on if the hit
                // is not on the current strand.
                if (bedItr->strand != strands[s])
                    continue;
                
                // new cluster, no overlap
                if ( ((bedItr->start - end) > _maxDistance) || (end < 0)) 
                {
                    cluster_id++;
                    end   = bedItr->end;
                }
                // same cluster, overlaps
                else {
                    if (bedItr->end > end) 
                        end = bedItr->end;
                }
                prev = *bedItr;
                _bed->reportBedTab(prev);
                printf("%d\n", cluster_id);
            }
        }
    }
}
