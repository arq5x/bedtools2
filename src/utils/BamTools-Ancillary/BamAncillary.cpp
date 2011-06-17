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
        if (bam.IsReverseStrand()) strand = "-";

        vector<CigarOp>::const_iterator cigItr = bam.CigarData.begin();
        vector<CigarOp>::const_iterator cigEnd = bam.CigarData.end();
        for ( ; cigItr != cigEnd; ++cigItr ) {
            if (cigItr->Type == 'M') {
                currPosition += cigItr->Length;
                blocks.push_back( BED(chrom, blockStart, currPosition, name, score, strand) );
                blockStart    = currPosition;
            }
            else if (cigItr->Type == 'D') {
                if (breakOnDeletionOps == false)
                    currPosition += cigItr->Length;
                else {
                    currPosition += cigItr->Length;
                    blockStart    = currPosition;
                }
            }
            else if (cigItr->Type == 'N') {
                currPosition += cigItr->Length;
                blockStart    = currPosition;            }
            else if (cigItr->Type == 'S' || cigItr->Type == 'H' || cigItr->Type == 'P' || cigItr->Type == 'I') {
                // do nothing
            }
            else {
                cerr << "Input error: invalid CIGAR type (" << cigItr->Type
                    << ") for: " << bam.Name << endl;
                exit(1);
            }
        }
    }
}
