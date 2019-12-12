/*****************************************************************************
  reldist.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "reldist.h"

/*
    Constructor
*/
RelativeDistance::RelativeDistance(string bedAFile, 
                                   string bedBFile,
                                   bool summary)
{
    _bedAFile  = bedAFile;
    _bedBFile  = bedBFile;
    _summary   = summary;
    _tot_queries = 0;
    CalculateRelativeDistance();
}


/*
    Destructor
*/
RelativeDistance::~RelativeDistance(void) {
}


void RelativeDistance::LoadMidpoints() {

    _bedB = new BedFile(_bedBFile);

    BED bed;    
    _bedB->Open();
    while (_bedB->GetNextBed(bed)) {
        CHRPOS midpoint = (int) (bed.end + bed.start) / 2;
        _db_midpoints[bed.chrom].push_back(midpoint);
    }
    
    map<string, vector<CHRPOS> >::const_iterator midItr = _db_midpoints.begin();
    map<string, vector<CHRPOS> >::const_iterator midEnd = _db_midpoints.end();
    for (; midItr != midEnd; ++midItr)
    {
        sort(_db_midpoints[midItr->first].begin(), 
             _db_midpoints[midItr->first].end());
    }
}


void RelativeDistance::ReportDistanceSummary()
{
    cout << "reldist\t"
         << "count\t"
         << "total\t"
         << "fraction\n";
    
    map<double, size_t>::const_iterator freqItr = _reldists.begin();
    map<double, size_t>::const_iterator freqEnd = _reldists.end();
    for (; freqItr != freqEnd; ++freqItr)
    {
        printf("%.2f\t%lu\t%lu\t%.3lf\n", 
               freqItr->first, 
               freqItr->second,
               _tot_queries,
               (double) freqItr->second / (double) _tot_queries);
    }
}


void RelativeDistance::UpdateDistanceSummary(double rel_dist)
{
    _tot_queries++;
    // round the relative distance to two decimal places.
    double rounded_rel_dist = floor(rel_dist * 100) / 100;
    _reldists[rounded_rel_dist]++;
}


void RelativeDistance::CalculateRelativeDistance()
{
    LoadMidpoints();
    
    vector<CHRPOS>::iterator low;
    CHRPOS low_idx, high_idx;
    double rel_dist;
    
    _bedA = new BedFile(_bedAFile);

    BED bed;    
    _bedA->Open();
    while (_bedA->GetNextBed(bed)) {
        
        if (_bedA->_status != BED_VALID)
            continue;

        if (_db_midpoints.count(bed.chrom) > 0)
        {
            vector<CHRPOS> *chrom_mids = &_db_midpoints[bed.chrom];
            // binary search the current query's midpoint among
            // the database midpoints
            CHRPOS midpoint = (CHRPOS) (bed.end + bed.start) / 2;
            low = lower_bound(chrom_mids->begin(), 
                              chrom_mids->end(), 
                              midpoint);
            
            // grab the indices for the database midpoints that are left and
            // right of the query's midpoint.
            if(low == chrom_mids->begin())
            {
                low_idx = low - chrom_mids->begin();
            }
            else 
            {
                low_idx = low - chrom_mids->begin() - 1;
            }
            high_idx = low_idx + 1;

            
            // make sure we don't run off the boundaries of the database's
            // midpoint vector
            if (low_idx != chrom_mids->size() - 1)
            {
                // grab the database midpoints that are left and right of
                // the query's midpoint.
                CHRPOS left = (*chrom_mids)[low_idx];
                CHRPOS right = (*chrom_mids)[high_idx];
                
                // ?
                if (left > midpoint)
                    continue;
                
                // calculate the relative distance between the query's midpoint
                // and the two nearest database midpoints.
                CHRPOS left_dist = abs(midpoint-left);
                CHRPOS right_dist = abs(midpoint-right); 
                if (min(left_dist, right_dist) == 0)
                {
                    rel_dist = 0.0;
                }
                else
                {
                    rel_dist = (double) min(left_dist, right_dist) / (double) (right-left);
                }
                if (!_summary)
                {
                    _bedA->reportBedTab(bed);
                    printf("%.3lf\n", rel_dist);
                }
                else { 
                    UpdateDistanceSummary(rel_dist);
                }
            }
        }
    }

    // report the "histogram" of distances.
    if (_summary)
        ReportDistanceSummary();
}


