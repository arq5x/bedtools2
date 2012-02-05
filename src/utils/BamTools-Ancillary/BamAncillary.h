/*****************************************************************************
  bamAncillary.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"
#include "lineFileUtilities.h"
#include "api/BamAlignment.h"

namespace BamTools {
    void getBamBlocks(const BamAlignment &bam, const RefVector &refs,
                        vector<BED> &blocks, bool includeDeletions = true);

    void MakeBedFromBam(const BamAlignment &bam, const string &chrom,
        const bedVector &blocks, BED &bed);
}
