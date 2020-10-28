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
    float l, r;

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
               l = _leftSlop;   
               _leftSlop  = _leftSlop * (float)bedEntry.size();
               r = _rightSlop;  
               _rightSlop = _rightSlop * (float)bedEntry.size();
               AddSlop(bedEntry);
               _rightSlop = r;
               _leftSlop = l;
           }
           _bed->reportBedNewLine(bedEntry);
       }
   }
   _bed->Close();
}
static inline CHRPOS _normalize_coord(CHRPOS size, CHRPOS pos) {
	if(pos < 0) return 0;
	if(size >= 0 && pos > size) return size;
	return pos;
}
void BedSlop::AddSlop(BED &bed) {

    // special handling if the BED entry is on the negative
    // strand and the user cares about strandedness.
    CHRPOS chromSize = (CHRPOS)_genome->getChromSize(bed.chrom);

	if(chromSize < 0) {
		cerr << "* Input error: Chromosome " << bed.chrom << " doesn't present in the .genome file. *" << endl;
		exit(1);
	}

	bool should_swap = _forceStrand && bed.strand == "-";
	CHRPOS left_slop = should_swap ? (CHRPOS)_rightSlop : (CHRPOS)_leftSlop;
	CHRPOS right_slop = should_swap ? (CHRPOS)_leftSlop : (CHRPOS)_rightSlop;

	bed.start -= left_slop;
	bed.end += right_slop;

	bed.start = _normalize_coord(chromSize, bed.start);
	bed.end = _normalize_coord(chromSize, bed.end);

    //checking edge case and adjusting
    if( bed.start == bed.end && bed.start == chromSize ){
      bed.start -= 1;
    }
    // Since bed files have 0 base system, the end is incremented by one
    else if ( bed.start == bed.end ){
      bed.end += 1;
    }
    else if(bed.start > bed.end){
      CHRPOS temp = bed.start;
      bed.start = bed.end;
      bed.end = temp;
    }
}
