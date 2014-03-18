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
                AddSlop(bedEntry);
            }
            else {
                _leftSlop  = _leftSlop  * (float)bedEntry.size();
                _rightSlop = _rightSlop * (float)bedEntry.size();
                AddSlop(bedEntry);
            }
            _bed->reportBedNewLine(bedEntry);
        }
    }
    _bed->Close();
}


void BedSlop::AddSlop(BED &bed) {

    // special handling if the BED entry is on the negative
    // strand and the user cares about strandedness.
    float chromSize = (float)_genome->getChromSize(bed.chrom);

    if ( (_forceStrand) && (bed.strand == "-") ) {
        // inspect the start
    	float newStart = (float)bed.start - _rightSlop;
    	bed.start = (newStart > 0 ) ? (CHRPOS)newStart : 0;

        // inspect the end
    	float newEnd = (float)bed.end + _leftSlop;
    	bed.end = (newEnd < chromSize ) ? (CHRPOS)newEnd : (CHRPOS)chromSize;
    }
    else {
        // inspect the start
    	float newStart = (float)bed.start - _leftSlop;
    	bed.start = (newStart > 0 ) ? (CHRPOS)newStart : 0;

        // inspect the end
    	float newEnd = (float)bed.end + _rightSlop;
    	bed.end = (newEnd < chromSize ) ? (CHRPOS)newEnd : (CHRPOS)chromSize;
    }
}


