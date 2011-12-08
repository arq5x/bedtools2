/*****************************************************************************
  slopBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "slopBed.h"


BedSlop::BedSlop(string &bedFile, string &genomeFile, bool forceStrand, 
                 float leftSlop, float rightSlop, bool fractional,
                 bool printHeader) {

    _bedFile     = bedFile;
    _genomeFile  = genomeFile;
    _forceStrand = forceStrand;
    _leftSlop    = leftSlop;
    _rightSlop   = rightSlop;
    _fractional  = fractional;
    _printHeader = printHeader;

    _bed    = new BedFile(bedFile);
    _genome = new GenomeFile(genomeFile);

    // get going, slop it up.
    SlopBed();
}


BedSlop::~BedSlop(void) {

}


void BedSlop::SlopBed() {

    BED bedEntry;     // used to store the current BED line from the BED file.

    _bed->Open();
    // report header first if asked.
    if (_printHeader == true) {
        _bed->PrintHeader();
    }        
    while (_bed->GetNextBed(bedEntry)) {    
        if (_bed->_status == BED_VALID) {
            if (_fractional == false) {
                AddSlop(bedEntry, (int) _leftSlop, (int) _rightSlop);
            }
            else {
                int leftSlop  = (int) (_leftSlop  * bedEntry.size());
                int rightSlop = (int) (_rightSlop * bedEntry.size());
                AddSlop(bedEntry, leftSlop, rightSlop);
            }
            _bed->reportBedNewLine(bedEntry);
        }
    }
    _bed->Close();
}


void BedSlop::AddSlop(BED &bed, int leftSlop, int rightSlop) {

    // special handling if the BED entry is on the negative
    // strand and the user cares about strandedness.
    CHRPOS chromSize = _genome->getChromSize(bed.chrom);

    if ( (_forceStrand) && (bed.strand == "-") ) {
        // inspect the start
        if ( (static_cast<int>(bed.start) - rightSlop) > 0 ) bed.start -= rightSlop;
        else bed.start = 0;

        // inspect the start
        if ( (static_cast<int>(bed.end) + leftSlop) <= static_cast<int>(chromSize)) bed.end += leftSlop;
        else bed.end = chromSize;
    }
    else {
        // inspect the start
        if ( (static_cast<int>(bed.start) - leftSlop) > 0) bed.start -= leftSlop;
        else bed.start = 0;

        // inspect the end
        if ( (static_cast<int>(bed.end) + rightSlop) <= static_cast<int>(chromSize)) bed.end += rightSlop;
        else bed.end = chromSize;
    }
}


