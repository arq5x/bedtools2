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
    float l = _leftSlop, r = _rightSlop;

    _bed->Open();
    // report header first if asked.
    if (_printHeader == true) {
        _bed->PrintHeader();
    }        
    while (_bed->GetNextBed(bedEntry)) {    
        if (_bed->_status == BED_VALID) {
            if (_fractional == false) {
                AddSlop(bedEntry);
            }
            else {
                _leftSlop  = l * bedEntry.size();
                _rightSlop = r * bedEntry.size();
                AddSlop(bedEntry);
            }
            _bed->reportBedNewLine(bedEntry);
        }
    }
	_leftSlop = l;
	_rightSlop = r;

    _bed->Close();
}


void BedSlop::AddSlop(BED &bed) {

    // special handling if the BED entry is on the negative
    // strand and the user cares about strandedness.
    CHRPOS chromSize = (CHRPOS)_genome->getChromSize(bed.chrom);
	
	CHRPOS leftShift, rightShift;

    if ( (_forceStrand) && (bed.strand == "-") ) {
        leftShift = (CHRPOS)_rightSlop;
        rightShift = (CHRPOS)_leftSlop;
    } else {
        leftShift = (CHRPOS)_leftSlop;
        rightShift = (CHRPOS)_rightSlop;
    }
    
    if (bed.start >= leftShift)
        bed.start -= leftShift;
    else
        bed.start = 0;
    
    if ((rightShift > chromSize) || (bed.end > chromSize - rightShift))
        bed.end = chromSize;
    else
        bed.end += rightShift;
}


