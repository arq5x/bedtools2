/*****************************************************************************
  multiBamCov.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "multiBamCov.h"
#include "api/BamMultiReader.h"


/*
    Constructor
*/
MultiCovBam::MultiCovBam(const vector<string> &bam_files, const string bed_file, 
                         int minQual, bool properOnly,
                         bool keepDuplicates, bool keepFailedQC,
                         bool obeySplits, bool sameStrand, 
                         bool diffStrand, float overlapFraction,
                         bool reciprocal)
:
_bam_files(bam_files),
_bed_file(bed_file),
_minQual(minQual),
_properOnly(properOnly),
_keepDuplicates(keepDuplicates),
_keepFailedQC(keepFailedQC),
_obeySplits(obeySplits),
_sameStrand(sameStrand),
_diffStrand(diffStrand),
_overlapFraction(overlapFraction),
_reciprocal(reciprocal)
{
	_bed = new BedFile(_bed_file);
    LoadBamFileMap();
}


/*
    Destructor
*/
MultiCovBam::~MultiCovBam(void) 
{}



bool MultiCovBam::FindBlockedOverlaps(const BED &a, const vector<BED> &a_blocks, 
                                      const BED &hit) {

    int a_footprint = GetTotalBlockLength(a_blocks);
    
    // 1. Break the raw hit into it's blocks
    //    and see of one of the hit blocks (inner loop)
    //    overlaps one of a's blocks (inner, inner loop)
    // 2. If so, mark the hit as valid and add it to the valid_set.
    //    Otherwise, the hit only overlapped the span of a and not
    //    the individual blocks.  Thus, it doesn't count.

    // break the hit into blocks
    bedVector hitBlocks;
    GetBedBlocks(hit, hitBlocks);
    int b_footprint = GetTotalBlockLength(hitBlocks);
    // test to see if there is a valid hit with one of the blocks
    bool valid_hit    = false;
    int tot_overlap = 0;
    bedVector::const_iterator hbItr = hitBlocks.begin();
    bedVector::const_iterator hbEnd = hitBlocks.end();
    for (; hbItr != hbEnd; ++hbItr) {
        // look for overlaps between this hit/block and each block in a
        bedVector::const_iterator a_blockItr = a_blocks.begin();
        bedVector::const_iterator a_blockEnd = a_blocks.end();
        for (; a_blockItr != a_blockEnd; ++a_blockItr) {
            int hs = max(a_blockItr->start, hbItr->start);
            int he = min(a_blockItr->end, hbItr->end);
            int overlap = he - hs;
            if (overlap > 0) {
                valid_hit = true;
                tot_overlap += overlap;
            }
        }
    }
    if (valid_hit)
    {
        // require sufficient overlap fraction (reciprocal or otherwise)
        // w.r.t to the "footprint" (i.e., the total length of each block)
        if ( ((float) tot_overlap / (float) a_footprint) > _overlapFraction) {
            if (_reciprocal && 
                ((float) tot_overlap / (float) b_footprint) > _overlapFraction) 
            {
                return true;
            }
            else if (!_reciprocal) {
                return true;
            }
        }
    }
    return false;
}


void MultiCovBam::CollectCoverage()
{
    BamMultiReader reader;
    
    if ( !reader.Open(_bam_files) )
    {
        cerr << "Could not open input BAM files." << endl;
        cerr << reader.GetErrorString() << endl;
        exit(1);
    }
    else
    {
        // attempt to find index files
        reader.LocateIndexes();

        // if index data available for all BAM files, we can use SetRegion
        if ( reader.HasIndexes() ) {
            BED bed;

            _bed->Open();
            // loop through each BED entry, jump to it, 
            // and collect coverage from each BAM
            while (_bed->GetNextBed(bed))
            {
                if (_bed->_status == BED_VALID)
                {
                    // initialize counts for each file to 0
                    vector<int> counts(_bam_files.size(), 0);
                    // get the BAM refId for this chrom.
                    int refId = reader.GetReferenceID(bed.chrom);
                    // set up a BamRegion to which to attempt to jump
                    BamRegion region(refId, (int)bed.start, refId, (int)bed.end);
                    
                    // everything checks out, just iterate through 
                    // specified region, counting alignments
                    if ( (refId != -1) && (reader.SetRegion(region)) ) {
                        BamAlignment al;
                        while ( reader.GetNextAlignment(al) )
                        {
                            string strand = "+";
                            if (al.IsReverseStrand() == true) strand = "-";
                            bool strands_are_same = (bed.strand == strand);

                            bool duplicate = al.IsDuplicate();
                            bool failedQC  = al.IsFailedQC();
                            if (_keepDuplicates) duplicate = false;
                            if (_keepFailedQC)    failedQC = false;

                            // filters
                            if (
                                (_properOnly && !al.IsProperPair()) ||
                                (_sameStrand && !strands_are_same)  ||
                                (_diffStrand && strands_are_same)   || 
                                (al.MapQuality < _minQual)          ||
                                (duplicate)                         ||
                                (failedQC)
                            )
                            {
                                continue;
                            }
                            

                            if (_obeySplits == false) {
                                // enforce fractional overlap
                                int al_end = al.GetEndPosition(false, false);
                                CHRPOS s = max((int)al.Position, (int) bed.start);
                                CHRPOS e = min(al_end, (int) bed.end);
                                CHRPOS aLength = (bed.end - bed.start);
                                CHRPOS bLength = (al_end - al.Position);
                                int overlapBases = (e - s);
                                float aOverlap = 
                                    ( (float) overlapBases / (float) aLength );
                                float bOverlap = 
                                    ( (float) overlapBases / (float) bLength );
                                
                                if ( aOverlap >= _overlapFraction) 
                                {
                                    if (!_reciprocal)
                                        counts[bamFileMap[al.Filename]]++;
                                    else if (bOverlap >= _overlapFraction)
                                        counts[bamFileMap[al.Filename]]++;
                                }
                            }
                            else {
                                // break alignment into discrete blocks,
                                bedVector bed_blocks, hits;
                                GetBamBlocks(al, bed.chrom, 
                                             bed_blocks, false, true);
                                // find the overlaps b/w the block in A & B
                                bool overlapsFound = FindBlockedOverlaps(bed, 
                                                                 bed_blocks, 
                                                                 bed);
                                if (overlapsFound == true)
                                    counts[bamFileMap[al.Filename]]++;
                            }
                        }
                    }
                    // report the cov at this interval for each file and reset
                    _bed->reportBedTab(bed);
                    ReportCounts(counts);
                }
            }
            _bed->Close();
        }
        else {
            cerr << "Could not find/load indexes." << endl;
            cerr << reader.GetErrorString() << endl;
            reader.Close();
            exit(1);
        }
    }
}


void MultiCovBam::LoadBamFileMap(void) 
{
    for (size_t i = 0; i < _bam_files.size(); ++i)
    {
        bamFileMap[_bam_files[i]] = i;
    }
}

void MultiCovBam::ReportCounts(const vector<int> &counts) 
{
    for (size_t i = 0; i < counts.size(); ++i)
    {
        if (i < counts.size() - 1)
            cout << counts[i] << "\t";
        else
            cout << counts[i];
    }
    cout << endl;
}
