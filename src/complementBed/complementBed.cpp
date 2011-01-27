/*****************************************************************************
  complementBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "complementBed.h"

BedComplement::BedComplement(string &bedFile, string &genomeFile) {

    _bedFile = bedFile;
    _genomeFile = genomeFile;

    _bed    = new BedFile(bedFile);
    _genome = new GenomeFile(genomeFile);

}


BedComplement::~BedComplement(void) {
}


//
// Merge overlapping BED entries into a single entry
//
void BedComplement::ComplementBed() {

    // load the "B" bed file into a map so
    // that we can easily compare "A" to it for overlaps
    _bed->loadBedFileIntoMapNoBin();

    // get a list of the chroms in the user's genome
    vector<string> chromList =  _genome->getChromList();

    // process each chrom in the genome
    for (size_t c = 0; c < chromList.size(); ++c) {
        string currChrom = chromList[c];

        // create a "bit vector" for the chrom
        CHRPOS currChromSize = _genome->getChromSize(currChrom);
        vector<bool> chromMasks(currChromSize, 0);

        // mask the chrom for every feature in the BED file
        bedVector::const_iterator bItr = _bed->bedMapNoBin[currChrom].begin();
        bedVector::const_iterator bEnd = _bed->bedMapNoBin[currChrom].end();
        for (; bItr != bEnd; ++bItr) {
            if (bItr->end > currChromSize) {
                cout << "Warninge: end of BED entry exceeds chromosome length. Please correct." << endl;
                _bed->reportBedNewLine(*bItr);
                exit(1);
            }

            // mask all of the positions spanned by this BED entry.
            for (CHRPOS b = bItr->start; b < bItr->end; b++)
                chromMasks[b] = 1;
        }

        // report the unmasked, that is, complemented parts of the chrom
        CHRPOS i = 0;
        CHRPOS start;
        while (i < chromMasks.size()) {
            if (chromMasks[i] == 0) {
                start = i;
                while ((chromMasks[i] == 0) && (i < chromMasks.size()))
                    i++;

                if (start > 0)
                    cout << currChrom << "\t" << start << "\t" << i << endl;
                else
                    cout << currChrom << "\t" << 0 << "\t" << i << endl;
            }
            i++;
        }
    }
}

