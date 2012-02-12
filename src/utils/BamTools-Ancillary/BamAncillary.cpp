/*****************************************************************************
  bamAncillary.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "BamAncillary.h"
using namespace std;

// 10   15   20      25    30               4000
// acccctttggacct---ataggga.................aaaa
// acccc---ggaccttttataggga.................aaaa
// 5M   3D 6M    2I 7M      20N             4M

namespace BamTools {
    void getBamBlocks(const BamAlignment &bam, const RefVector &refs,
                      vector<BED> &blocks, bool breakOnDeletionOps) {

        CHRPOS currPosition = bam.Position;
        CHRPOS blockStart   = bam.Position;
        string chrom        = refs.at(bam.RefID).RefName;
        string name         = bam.Name;
        string strand       = "+";
        string score        = ToString(bam.MapQuality);
        char  prevOp        = '\0';
        if (bam.IsReverseStrand()) strand = "-";
        bool blocksFound = false;

        vector<CigarOp>::const_iterator cigItr = bam.CigarData.begin();
        vector<CigarOp>::const_iterator cigEnd = bam.CigarData.end();
        for ( ; cigItr != cigEnd; ++cigItr ) {
            if (cigItr->Type == 'M') {
                currPosition += cigItr->Length;
                // we only want to create a new block if the current M op
                // was preceded by an N op or a D op (and we are breaking on D ops)
                if ((prevOp == 'D' && breakOnDeletionOps == true) || (prevOp == 'N')) {
                    blocks.push_back( BED(chrom, blockStart, currPosition, name, score, strand) );
                    blockStart = currPosition;
                }
            }
            else if (cigItr->Type == 'D') {
                if (breakOnDeletionOps == false)
                    currPosition += cigItr->Length;
                else {
                    blocksFound = true;
                    currPosition += cigItr->Length;
                    blockStart    = currPosition;
                }
            }
            else if (cigItr->Type == 'N') {
                blocks.push_back( BED(chrom, blockStart, currPosition, name, score, strand) );
                blocksFound = true;
                currPosition += cigItr->Length;
                blockStart    = currPosition;
            }
            else if (cigItr->Type == 'S' || cigItr->Type == 'H' || cigItr->Type == 'P' || cigItr->Type == 'I') {
                // do nothing
            }
            else {
                cerr << "Input error: invalid CIGAR type (" << cigItr->Type
                    << ") for: " << bam.Name << endl;
                exit(1);
            }
            prevOp = cigItr->Type;
        }
        // if there were no splits, we just create a block representing the contiguous alignment.
        if (blocksFound == false) {
            blocks.push_back( BED(chrom, bam.Position, currPosition, name, score, strand) );
        }
    }
    
    void MakeBedFromBam(const BamAlignment &bam, const string &chrom,
                        const bedVector &blocks, BED &bed) 
    {
        bed.chrom = chrom;
        bed.start = bam.Position;
        bed.end   = bam.GetEndPosition(false, false);
        // build the name field from the BAM alignment.
        bed.name = bam.Name;
        if (bam.IsFirstMate()) bed.name += "/1";
        else if (bam.IsSecondMate()) bed.name += "/2";
        bed.score  = ToString(bam.MapQuality);
        bed.strand = "+";
        if (bam.IsReverseStrand()) bed.strand = "-";
        
        bed.fields.push_back(bed.chrom);
        bed.fields.push_back(ToString(bed.start));
        bed.fields.push_back(ToString(bed.end));
        bed.fields.push_back(bed.name);
        bed.fields.push_back(bed.score);
        bed.fields.push_back(bed.strand);
        bed.fields.push_back(ToString(bam.Position));
        bed.fields.push_back(ToString(bed.end));
        bed.fields.push_back("0,0,0");
        bed.fields.push_back(ToString(blocks.size()));
        
        ostringstream block_lengths, block_starts;
        for (size_t i = 0; i < blocks.size(); ++i) {
            block_lengths << blocks[i].end   - blocks[i].start << ",";
            block_starts  << blocks[i].start - bam.Position << ",";
        }
        bed.fields.push_back(block_lengths.str());
        bed.fields.push_back(block_starts.str());
        
        // identify which indices relate to the "other" 
        // (i,e. non-BED6) fields
        for (size_t i = 6; i <= 11; ++i)
            bed.other_idxs.push_back(i);
    }
}
