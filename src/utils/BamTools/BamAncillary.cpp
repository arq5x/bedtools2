/*****************************************************************************
  bamAncillary.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "BamAncillary.h"
using namespace std;

// 10   15   20      25    30               4000
// acccctttggacct---ataggga.................aaaa
// acccc---ggaccttttataggga.................aaaa
// 5M   3D 6M    2I 7M      20N             4M
  
namespace BamTools {
    void getBamBlocks(const BamAlignment &bam, const RefVector &refs, 
                      vector<BED> &blocks, bool includeDeletions) {
    
    	CHRPOS currPosition = bam.Position;
        CHRPOS blockStart   = bam.Position;
        string chrom        = refs.at(bam.RefID).RefName;

    	vector<CigarOp>::const_iterator cigItr = bam.CigarData.begin();
    	vector<CigarOp>::const_iterator cigEnd = bam.CigarData.end();
        for ( ; cigItr != cigEnd; ++cigItr ) {
    		switch (cigItr->Type) {
    		case 'M':
                currPosition += cigItr->Length;
    			blocks.push_back( BED(chrom, blockStart, currPosition) );
    			break;

    		case 'S':
    		case 'P':
    		case 'H':
    		case 'I':
    			// Insertion relative to the reference genome.
    			// Don't advance the current position, since no new nucleotides are covered.
    			break;

    		case 'D':
    		    if (includeDeletions == true)
    		        currPosition += cigItr->Length;
    		    else {
                    blocks.push_back( BED(chrom, blockStart, currPosition) );
                    currPosition += cigItr->Length;
                    blockStart    = currPosition;		        
    		    }
    		case 'N':
                currPosition += cigItr->Length;
                blockStart    = currPosition;            
    			break;

    		default:
    			cerr << "Input error: invalid CIGAR type (" << cigItr->Type
    				<< ") for: " << bam.Name << endl;
    			exit(1);
    		}
    	}
    }
}