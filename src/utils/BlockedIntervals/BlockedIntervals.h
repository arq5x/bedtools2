/*****************************************************************************
  BlockedIntervals.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include <vector>
#include "bedFile.h"
#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

/*
    Using the CIGAR string, break a BAM alignment
    into discrete alignment blocks.
*/
void GetBamBlocks(const BamAlignment &bam,
                      const string &chrom,
                      bedVector &bedBlocks,
                      bool breakOnDeletionOps,
                      bool breakOnSkipOps);

/* break a BED12 record into discrete BED6 blocks. */
void GetBedBlocks(const BED &bed, bedVector &bedBlocks);

/* compute the total forprint of a set of BED blocks */
int GetTotalBlockLength(const bedVector &bedBlocks);
