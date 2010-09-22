/*****************************************************************************
  slopBed.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "slopBed.h"


BedSlop::BedSlop(string &bedFile, string &genomeFile, bool &forceStrand, int &leftSlop, int &rightSlop) {

	_bedFile = bedFile;
	_genomeFile = genomeFile;
	_forceStrand = forceStrand;
	
	_leftSlop = leftSlop;
	_rightSlop = rightSlop;
	
	_bed    = new BedFile(bedFile);
	_genome = new GenomeFile(genomeFile);
	
	SlopBed();
}


BedSlop::~BedSlop(void) {

}


void BedSlop::SlopBed() {
	
	int lineNum = 0;
	BED bedEntry, nullBed;     // used to store the current BED line from the BED file.
	BedLineStatus bedStatus;
	
	_bed->Open();
	bedStatus = _bed->GetNextBed(bedEntry, lineNum);
	while (bedStatus != BED_INVALID) {
		if (bedStatus == BED_VALID) {
			AddSlop(bedEntry);
			_bed->reportBedNewLine(bedEntry);
			bedEntry = nullBed;	
		}
		bedStatus = _bed->GetNextBed(bedEntry, lineNum);				
	}
	_bed->Close();
}


void BedSlop::AddSlop(BED &bed) {

	// special handling if the BED entry is on the negative
	// strand and the user cares about strandedness.
	CHRPOS chromSize = _genome->getChromSize(bed.chrom);
	
	if ( (_forceStrand) && (bed.strand == "-") ) {
		// inspect the start
		if ( (static_cast<int>(bed.start) - _rightSlop) > 0 ) bed.start -= _rightSlop;
		else bed.start = 0;

		// inspect the start		
		if ( (static_cast<int>(bed.end) + _leftSlop) <= static_cast<int>(chromSize)) bed.end += _leftSlop;
		else bed.end = chromSize;
	}
	else {		
		// inspect the start
		if ( (static_cast<int>(bed.start) - _leftSlop) > 0) bed.start -= _leftSlop;
		else bed.start = 0;
		
		// inspect the end
		if ( (static_cast<int>(bed.end) + _rightSlop) <= static_cast<int>(chromSize)) bed.end += _rightSlop;
		else bed.end = chromSize;
	}
}


